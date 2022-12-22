# Differential Equations

## Euler Method
```@docs
euler(F::Function,a::Number,y_0::Number,b::Number,nh::Number)
```

## Modified Euler Method
```@docs
modified_euler(F::Function,a::Number,y_0::Number,b::Number,nh::Number)
```

## Runge Kutta of Order 2 Method
```@docs
runge_kutta_order_2(F::Function,a::Number,y_0::Number,b::Number,nh::Number)
```

## Heun Method
```@docs
heun(F::Function,a::Number,y_0::Number,b::Number,nh::Number)
```