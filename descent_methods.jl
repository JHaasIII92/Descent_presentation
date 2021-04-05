function GD_ellispe(∇f::Function, b, x_k, y_k, ϵ)
    k = 1
    x_vals = [x_k]
    y_vals = [y_k]
    while norm(∇f(x_k,y_k)) > ϵ
        x_k = b*((b-1)/(b+1))^k
        y_k = ((1-b)/(b+1))^k
        k+=1
        
        push!(x_vals, x_k)
        push!(y_vals, y_k)
    end
    return x_vals, y_vals
end

function Quad_Backtracking(f::Function, ∇f::Function, P, x_k, y_k, α, β, ϵ)

    x_vals = []
    y_vals = []
    
    ∇fx = ∇f(x_k, y_k)
    
    for i in 1:200
        if norm(∇fx) < ϵ
            break
        end
        ∇fx = ∇f(x_k, y_k)
        fx = f(x_k, y_k)
        Δxk = P\∇fx
        ∇fxT_Δxk = ∇fx'Δxk
        t = 1
        while f(x_k + t*Δxk[1], y_k + t*Δxk[2]) > fx + α*t*∇fxT_Δxk
            t = β * t
        end
        x_k = x_k + t*Δxk[1]
        y_k = y_k + t*Δxk[2]
        
        push!(x_vals, x_k)
        push!(y_vals, y_k)
    end
    return x_vals, y_vals
end 

function updating_Hessien_Descent(f::Function, ∇f::Function, H::Function, x_k, y_k, l, α, β, ϵ)
    
    
    x_vals = []
    y_vals = []
    
    P = -1 * H(x_k, y_k)
    for i in 1:200
        ∇fx = ∇f(x_k, y_k)
        if norm(∇fx) < ϵ
            break
        end
        
        if i % l == 0
            P = -1 * H(x_k, y_k)
        end
        
        fx = f(x_k, y_k)
        Δxk = P\∇fx
        ∇fxT_Δxk = ∇fx'Δxk
        while f(x_k + t*Δxk[1], y_k + t*Δxk[2]) > fx + α*t*∇fxT_Δxk
            t = β * t
        end
        
        x_k = x_k + t*Δxk[1]
        y_k = y_k + t*Δxk[2]
        
        push!(x_vals, x_k)
        push!(y_vals, y_k)
    end
    return x_vals, y_vals
end 