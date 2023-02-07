using GLMakie
using ProgressBars

function calcGrad(vx, vy)
    _size = size(vx)
    grad_x = zeros(Float32, _size)
    grad_y = zeros(Float32, _size)

    for x in 2:size(grad_x,2) - 1
        for y in 2:size(grad_x,1) - 1
            grad_x[y,x] = ((vx[y-1, x-1] + vx[y+1, x+1] + vx[y+1, x-1] + vx[y-1, 		x+1]) * 1 +
                (vx[y, x-1] + vx[y, x+1] - vx[y-1, x] - vx[y+1, x]) * 2 +
                vx[y,x] * -4) / 16

            grad_y[y,x] = ((vy[y-1, x-1] + vy[y+1, x+1] - vy[y+1, x-1] - vy[y-1, 		x+1]) * -1 +
                (vy[y, x-1] + vy[y, x+1] - vy[y-1, x] - vy[y+1, x]) * -2 +
                vy[y,x] * -4) / 16
        end
    end
    return (grad_x, grad_y)
end

mutable struct Fluid
    width::Int
    height::Int

    image::AbstractArray{UInt8}

    pressure::AbstractArray{Float32}
    density::AbstractArray{Float32}
    velocity_x::AbstractArray{Float32}
    velocity_y::AbstractArray{Float32}

    workspace_1::AbstractArray{Float32}
    workspace_2::AbstractArray{Float32}

    function Fluid(width, height)
        image     = zeros(UInt8, height, width)

        pressure    = zeros(Float32, height, width)
        density     = zeros(Float32, height, width)
        velocity_x  = zeros(Float32, height, width)
        velocity_y  = zeros(Float32, height, width)
        workspace_1 = zeros(Float32, height, width)
        workspace_2 = zeros(Float32, height, width)

        return new(width, height, image, pressure, density,
                    velocity_x, velocity_y, workspace_1, workspace_2)
    end
end

function linInterp(quantity, x, y)

    width = size(quantity,  2)
    height = size(quantity, 1)

    x1 = min(max(Int(floor(x)), 1), width)
    x2 = min(max( Int(ceil(x)), 1), width)

    y1 = min(max(Int(floor(y)), 1), height)
    y2 = min(max( Int(ceil(y)), 1), height)

    f11 = quantity[y1, x1]
    f12 = quantity[y1, x2]
    f21 = quantity[y2, x1]
    f22 = quantity[y2, x2]

    x = modf(x)[1]
    y = modf(y)[1]

    fxy1 = (1.0 - x) * f11 + x * f12
    fxy2 = (1.0 - x) * f21 + x * f22

    fxy = (1.0 - y) * fxy1 + y * fxy2

    # @info "lininterp" x y
    return fxy
end

function addInflow(fluid::Fluid, x, y, w, h, d, u, v)
    @assert x <= 1.0 && x >= 0.0 "x must be on [0.0, 1.0]"
    @assert y <= 1.0 && y >= 0.0 "y must be on [0.0, 1.0]"
    @assert w <= 1.0 && w > 0.0  "w must be on (0.0, 1.0]"
    @assert h <= 1.0 && h > 0.0  "h must be on (0.0, 1.0]"
    @assert x + w <= 1.0 "x + w can't exceed 1.0"
    @assert y + h <= 1.0 "y + h can't exceed 1.0"

    start_x = Int(round(x * (fluid.width - 1) + 1))
    start_y = Int(round(y * (fluid.height - 1) + 1))

    end_x = Int(round(start_x + w * (fluid.width - 1) ))
    end_y = Int(round(start_y + h * (fluid.height - 1) ))

    fluid.density[start_y:end_y, start_x:end_x] .= d

    fluid.velocity_x[start_y:end_y, start_x:end_x] .= u
    fluid.velocity_y[start_y:end_y, start_x:end_x] .= v

    return nothing
end

# Convert density to grayscale output
function densityImage!(image, density)
    image .= UInt8.(floor.(density .* 255.0))
    image .= max.(min.(image, 255), 0)
    return nothing
end

function advectDensity(fluid::Fluid, timestep)
    advectDensity(fluid.velocity_x, fluid.velocity_y, fluid.density,
                  fluid.workspace_1, timestep)
    return nothing
end

function advectDensity(v_x, v_y, density, density_new, timestep)
    height, width = size(density)

    for x in 1:width
        for y in 1:height
            new_x = x - v_x[y, x] * timestep
            new_y = y - v_y[y, x] * timestep
            density_new[y, x] = linInterp(density, new_x, new_y)
        end
    end

    density .= density_new

    return nothing
end

function advectVelocity(v_x, v_y, w_1, w_2, timestep)
    width, height = size(v_x)

    for x in 1:width
        for y in 1:height
            old_x = x - v_x[y, x] * timestep
            old_y = y - v_y[y, x] * timestep
            w_1[y, x] = linInterp(v_x, old_x, old_y)
            w_2[y, x] = linInterp(v_y, old_x, old_y)
        end
    end

    # Zero velocity at boundaries
    w_1[:,1]   .= 0
    w_1[:,end] .= 0

    w_2[1,:]   .= 0
    w_2[end,:] .= 0

    # Update velocity values
    v_x .= w_1
    v_y .= w_2

    return nothing
end

# Advect velocity values
function advectVelocity(fluid::Fluid, timestep)
    advectVelocity(fluid.velocity_x, fluid.velocity_y,
                    fluid.workspace_1, fluid.workspace_2,
                    timestep)

    return nothing
end

# Calculate the divergence of the velocity field
function velocityDivergenceDiff(fluid::Fluid)
    div_x = diff(fluid.velocity_x, dims=2)
    div_y = diff(fluid.velocity_y, dims=1)

    fluid.workspace_1 .= 0.0
    fluid.workspace_1[:,1:end-1] .+= div_x
    fluid.workspace_1[1:end-1, :] .+= div_y
    return nothing
end

# Solve for pressure distribution to eliminate divergence
function pressureSolve(pressure, divergence, workspace, max_iters)
    _size = size(pressure)

    # Jacobi iteration
    for i in 1:max_iters - 1
        for x in 2:_size[2] - 1
            for y in 2:_size[1] - 1
                workspace[y, x] =
                    (pressure[y, x+1] + pressure[y, x-1] + pressure[y+1, x] + pressure[y-1, x] -
                    divergence[y,x]) / 4
            end
        end
        # Swap pressure and workspace
        pressure, workspace = workspace, pressure
    end
    return nothing
end

# Iteratively solve for pressure
function pressureSolve(fluid::Fluid, num_iters)
    pressureSolve(fluid.pressure, fluid.workspace_1, fluid.workspace_2, num_iters)
    return nothing
end

# Diffuse density and velocity TODO
function diffuse(v_x, v_y, w_1, w_2, timestep, max_iters)
    # Initial guess of identity
    # pressure = deepcopy()

    # Jacobi iteration
    # for i in 1:num_iters - 1
        # for x in 2:fluid.width - 1
    #        for y in 2:fluid.height - 1
                # d_X = (d0_X + diffusionFactor * deltaTime * (d_01 + d_02+ d_03 + d_04)) /
                    # (1 + 4 * diffusionFactor * deltaTime)
    #        end
    #    end

    # end
end

function diffuse(fluid::Fluid, timestep, max_iters)
    diffuse(fluid.velocity_x, fluid.velocity_y,
            fluid.workspace_1, fluid.workspace_2, timestep, max_iters)
    return nothing
end

function calcDivergence(vx, vy, div)
    _size = size(vx)

    for x in 2:_size[2] - 1
        for y in 2:_size[1] - 1
            div[y,x] = (vy[y+1, x] - vy[y,x]) - (vy[y-1, x] - vy[y,x]) + (vx[y, x+1] - vx[y,x]) - (vx[y, x-1] - vx[y,x])
        end
    end
    return nothing
end


# Correct the cell velocities to eliminate divergence
function correctVel(v_x, v_y, p, v_x_new, v_y_new, scale)
    _size = size(v_x)

    for x in 2:_size[2] - 1
        for y in 2:_size[1] - 1
            v_x_new[y,x] = v_x[y,x] - (p[y,x] - p[y, x-1]) * scale
            v_y_new[y,x] = v_y[y,x] - (p[y,x] - p[y-1, x]) * scale
        end
    end

    v_x_new[:,1]   .= 0.0
    v_x_new[:,end] .= 0.0
    v_y_new[1,:]   .= 0.0
    v_y_new[end,:] .= 0.0

    v_x .= v_x_new
    v_y .= v_y_new

    return nothing
end

function correctVel(fluid::Fluid, scale=0.1)
    correctVel(fluid.velocity_x, fluid.velocity_y, fluid.pressure,
               fluid.workspace_1, fluid.workspace_2, scale)
    return nothing
end

function applyPressure(fluid::Fluid, timestep, max_iters=1000)
    # Calculate divergence of current velocity field
    calcDivergence(fluid.velocity_x, fluid.velocity_y, fluid.workspace_1),

    # Solve for Pressure
    pressureSolve(fluid, max_iters),

    # Correct velocity filed to eliminate divergence
    correctVel(fluid.velocity_x, fluid.velocity_y, fluid.pressure,
               fluid.workspace_1, fluid.workspace_2, 1)
    return nothing
end

function step(fluid::Fluid, timestep, use_gpu::Bool=false)
    # Advect
    advectDensity(fluid, timestep)
    advectVelocity(fluid, timestep)

    # Pressure Solve
    applyPressure(fluid, timestep)

    # TODO
    # diffuse(fluid, timestep)
    return nothing
end

function simulate(width, height, sim_duration, timestep)
    fluid = Fluid(width, height)
    timestamps = collect(0.0:timestep:sim_duration)
    frames = Array{UInt8}(undef, (height, width, size(timestamps)[1]))

    x = collect(1:width)
    y = collect(1:height)

    coordinates_x = [x for x in x, y in y]
    coordinates_y = [y for x in x, y in y]

    force = 300.0 .* (
        20 .* exp.( - 1.0 / (2 * 0.005) * ((coordinates_x ./ width .- 0.45).^2 + (coordinates_y ./ height .- 0.30).^2))
        -
        exp.( - 1.0 / (2 * 0.005) * ((coordinates_x ./ width .- 0.55).^2 + (coordinates_y ./ height .- 0.70).^2)))

    # Gaussian density distribution to start
    fluid.density .= exp.( - 1.0 / (2 * 0.1) * ((coordinates_x ./ width .- 0.5).^2 + (coordinates_y ./ height .- 0.5).^2))

    # Run through simulation time
    progress_frames = ProgressBar(1:size(timestamps)[1])
    for i in progress_frames
        # Add inflow
        addInflow(fluid, 0.4, 0.8999999, 0.05, 0.1, 0.3, 0.0, -4000.0)
        addInflow(fluid, 0.4, 0.01, 0.05, 0.1, 0.3, 0.0, 4000.0)

        # pre_factor = max(2 - cur_time, 0.0)
        # fluid.velocity_x += timestep * pre_factor * force

        # Fluid Update
        step(fluid, timestep)

        # Update frame
        densityImage!(fluid.image, fluid.density)
        frames[:,:, i] .= fluid.image
   end

    return frames
end

function main()
    @info "Simulating..."
    frames = simulate(2000, 2000, 4.0, 0.01);

    @info "Displaying..."
    iter = ProgressBar(1:size(frames)[3])

    figure, ax, plot = heatmap(frames[:,:,1],
                                label=false,
                                xticks=false,
                                yticks=false,
                                colormap = :heat)

    hidedecorations!(ax)
    hidespines!(ax)

    # TODO: Wrap into simulation and make displayed array configurable
    # Generate output visualization media
    record(figure, "sim_4s_2k.mp4", iter) do i
        plot[1] = frames[:,:,i]
    end
end

isinteractive() || main()
