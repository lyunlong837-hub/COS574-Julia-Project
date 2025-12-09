# nonlinear_calculator.jl
#
# Implementation and Evaluation of Nonlinear Equation Solvers by Julia
#
# Methods:
#   - Bisection
#   - Newton (with numerical derivative if not provided)
#   - Secant
#   - Broyden (1D quasi-Newton)


# Utility: parse function
function make_function(expr::String)
    ex = Meta.parse(expr)
    f(x) = eval(:(let x = $x; $ex end))
    return f
end


# Bisection Method
function bisection(f, a::Float64, b::Float64; tol=1e-8, maxiter=100)
    fa = f(a)
    fb = f(b)
    if fa * fb > 0
        error("Bisection method requires f(a) * f(b) ≤ 0. Current: f(a) = $fa, f(b) = $fb")
    end

    history = Float64[]
    left, right = a, b
    mid = (left + right) / 2
    fmid = f(mid)

    for k in 1:maxiter
        mid = (left + right) / 2
        fmid = f(mid)
        push!(history, mid)

        if abs(fmid) < tol || (right - left)/2 < tol
            return mid, (converged=true, iterations=k, history=history)
        end

        if fa * fmid ≤ 0
            right = mid
            fb = fmid
        else
            left = mid
            fa = fmid
        end
    end

    return mid, (converged=false, iterations=maxiter, history=history)
end


# Newton Method
function newton(f, x0::Float64; df=nothing, tol=1e-8, maxiter=100)
    history = Float64[]
    x = x0

    # Numerical derivative if not provided
    function num_derivative(x)
        h = 1e-6 * max(1.0, abs(x))
        return (f(x + h) - f(x - h)) / (2h)
    end

    for k in 1:maxiter
        fx = f(x)
        push!(history, x)

        if abs(fx) < tol
            return x, (converged=true, iterations=k, history=history)
        end

        dfx = df === nothing ? num_derivative(x) : df(x)
        if dfx == 0
            return x, (converged=false, iterations=k, history=history)
        end

        x_new = x - fx / dfx

        if abs(x_new - x) < tol
            push!(history, x_new)
            return x_new, (converged=true, iterations=k, history=history)
        end

        x = x_new
    end

    return x, (converged=false, iterations=maxiter, history=history)
end


# Secant Method
function secant(f, x0::Float64, x1::Float64; tol=1e-8, maxiter=100)
    history = Float64[]
    f0 = f(x0)
    f1 = f(x1)
    push!(history, x0)
    push!(history, x1)

    for k in 1:maxiter
        if f1 == f0
            return x1, (converged=false, iterations=k, history=history)
        end

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        push!(history, x2)

        if abs(x2 - x1) < tol || abs(f(x2)) < tol
            return x2, (converged=true, iterations=k, history=history)
        end

        x0, f0 = x1, f1
        x1, f1 = x2, f(x2)
    end

    return x1, (converged=false, iterations=maxiter, history=history)
end

# Broyden Method (1D)
function broyden(f, x0::Float64, x1::Float64; tol=1e-8, maxiter=100)
    history = Float64[]
    f0 = f(x0)
    f1 = f(x1)
    push!(history, x0)
    push!(history, x1)

    # Initial approximate derivative
    if x1 == x0
        B = 1.0
    else
        B = (f1 - f0) / (x1 - x0)
    end

    x = x1
    fx = f1

    for k in 1:maxiter
        if B == 0
            return x, (converged=false, iterations=k, history=history)
        end

        s = -fx / B
        x_new = x + s
        fx_new = f(x_new)
        push!(history, x_new)

        if abs(fx_new) < tol || abs(x_new - x) < tol
            return x_new, (converged=true, iterations=k, history=history)
        end

        # Broyden update for approximate derivative B:
        # B_{k+1} = B_k + (y - B_k s) / s, where y = f(x_{k+1}) - f(x_k)
        y = fx_new - fx
        if s != 0
            B = B + (y - B * s) / s
        end

        x, fx = x_new, fx_new
    end

    return x, (converged=false, iterations=maxiter, history=history)
end






function prompt_float(msg::String)
    print(msg)
    return parse(Float64, readline())
end

function main()
    println("==============================================")
    println(" Nonlinear Equation Solver Calculator (Julia)")
    println(" Methods: Bisection, Newton, Secant, Broyden")
    println("==============================================")
    println()
    println("Please input f(x) as a Julia expression, e.g.")
    println("    x^3 - x - 2")
    println("    sin(x) - 0.5")
    println()

    print("f(x) = ")
    expr = readline()
    f = make_function(expr)

    println()
    println("Choose method")
    println("  1) Bisection")
    println("  2) Newton")
    println("  3) Secant")
    println("  4) Broyden")
    print("Your choice = ")
    method = parse(Int, readline())

    tol = prompt_float("Tolerance (default 1e-8, just press Enter to use default) or input a value: ")
    
    maxiter = Int(prompt_float("Max iterations (e.g. 100): "))

    println()

    if method == 1
        println("Bisection method requires an interval [a, b] with f(a)*f(b) ≤ 0.")
        a = prompt_float("a = ")
        b = prompt_float("b = ")

        root, info = bisection(f, a, b; tol=tol, maxiter=maxiter)
        println("\n===== Result (Bisection) =====")
        println("Root       ≈ $root")
        println("f(root)    = $(f(root))")
        println("Iterations = $(info.iterations)")
        println("Converged  = $(info.converged)")

    elseif method == 2
        println("Newton method uses one initial guess x0.")
        x0 = prompt_float("x0 = ")

        root, info = newton(f, x0; tol=tol, maxiter=maxiter)
        println("\n===== Result (Newton) =====")
        println("Root       ≈ $root")
        println("f(root)    = $(f(root))")
        println("Iterations = $(info.iterations)")
        println("Converged  = $(info.converged)")

    elseif method == 3
        println("Secant method uses two initial guesses x0, x1.")
        x0 = prompt_float("x0 = ")
        x1 = prompt_float("x1 = ")

        root, info = secant(f, x0, x1; tol=tol, maxiter=maxiter)
        println("\n===== Result (Secant) =====")
        println("Root       ≈ $root")
        println("f(root)    = $(f(root))")
        println("Iterations = $(info.iterations)")
        println("Converged  = $(info.converged)")

    elseif method == 4
        println("Broyden method (1D) uses two initial guesses x0, x1.")
        x0 = prompt_float("x0 = ")
        x1 = prompt_float("x1 = ")

        root, info = broyden(f, x0, x1; tol=tol, maxiter=maxiter)
        println("\n===== Result (Broyden) =====")
        println("Root       ≈ $root")
        println("f(root)    = $(f(root))")
        println("Iterations = $(info.iterations)")
        println("Converged  = $(info.converged)")

    else
        println("Invalid method choice.")
    end
end

# Only run main() when this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
