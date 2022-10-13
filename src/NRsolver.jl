function NRsolver(f::Function, df::Function, x0::Array{Float64}, reltol::Float64=1e-5, abstol::Float64=1e-5, n_max::Integer=20, jacobi="update")
    """
    solver of Newton-Raphson method 
    to solve funtion of 'f(x)=0' using iteration method:
    x_new = x_old - df(x_old)^(-1)*f(x_old)
    """

    x_old = x0
    r_old = f(x_old)
    converge = false
    for i = 1:n_max
        if jacobi == "update" && i <= 3 # only the beginning 3 steps update the jacobi matrix
            df_inv = inv(df(x_old))
        end
        x_new = x_old - df_inv * r_old
        r_new = f(x_new)
        if norm(x_new - x_old) < reltol || norm(r_new) < abstol
            converge = true
            return x_new
        end
    end
    if i == n_max && converge == false
        error("the maximum number of interation is reached and the result is not converged")
    end
end