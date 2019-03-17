using ODE;

C_m = 1;
#V_m = -58.7085;

I_app = 0;
I_stim = 0;

g_Na = 40;#120;
g_K = 35;#36;
g_L = 0.3;

E_Na = 55;
E_K = -77;
E_L = -65;#-54.4;
AverageFrequency = zeros(Float64, 0)
I_appArray = zeros(Float64, 0)

#r0 = [-58.7085; 0.0953; 0.000913; 0.3662];
#r0 = [14.8409; 0.9174;  0.0140; 0.0539]
r0 = [-20; 0.5;  0.5; 0.05]
const dt = 0.000025;
const tf  = 20.0;

#t = 0:dt:tf;

using Findpeaks

function frequency(V,t)
    Peaks = findpeaks(V, minHeight= -10.0)
    if length(Peaks) == 0 || length(Peaks) == 1
        return 0       #Frequency is null
    end

    Half = convert(Int, floor(length(Peaks)/2))
    sort!(Peaks)
    Peaks = Peaks[Half:end]
    Frequency = Array{Float64}(undef, length(Peaks)-1)

    for i = 1:length(Peaks)-1
        Frequency[i] = 1/(t[Peaks[i + 1]] - t[Peaks[i]])
    end

    return Frequency
end;

function f(t, r)
    (V, m, n, h) = r;

    dV_dt = 1000 * ((g_Na * m^3 * h * (E_Na - V)) + (g_K * n * (E_K - V)) + (g_L * (E_L - V)) + I_app);
    dm_dt = 1000 * (0.182 * (V + 35) / (1 - exp(-(V + 35) / 9)) * (1 - m) - (-0.124 * (V + 35) / (1 - exp((V + 35) / 9))) * m);
    dn_dt = 1000 * (0.02 * (V - 25) / ( 1 - exp(-(V - 25) / 9)) * (1 - n) - (-0.002 * (V - 25) / (1 - exp((V - 25) / 9))) * n);
    dh_dt = 1000 * (0.25 * exp(-(V + 90) / 12) * (1 - h) - 0.25 * exp((V + 62) / 6) / exp((V + 90) / 12) * h);

    [dV_dt; dm_dt; dn_dt; dh_dt];
end;

#for I_app in 0:0.1:1.8
#    t = 0:dt:tf;
#    println(r0)
#    println("hyu")
#    (t, r) = ode45(f, r0, t)
#    V = map(v -> v[1], r);
#    m = map(v -> v[2], r);
#    n = map(v -> v[3], r);
#    h = map(v -> v[4], r);
#    freq = frequency(V, t)

#    length(freq) == 0 ? push!(AverageFrequency, 0) : push!(AverageFrequency, sum(freq)/length(freq))
    #r0 = [V[length(V)], m[length(m)], n[length(n)], h[length(h)]]
#    empty!(r0)
#    push!(r0, V[length(V)], m[length(m)], n[length(n)], h[length(h)])
#    push!(I_appArray, I_app)
#end;

while I_app <= 1.9
    t = 0:dt:tf;
    println(r0, " r0")
    (t, r) = ode45(f, r0, t)
    V = map(v -> v[1], r);
    m = map(v -> v[2], r);
    n = map(v -> v[3], r);
    h = map(v -> v[4], r);
    freq = frequency(V, t)

    length(freq) == 0 ? push!(AverageFrequency, 0) : push!(AverageFrequency, sum(freq)/length(freq))
    println(AverageFrequency, " AverageFrequency")
    #r0 = [V[length(V)], m[length(m)], n[length(n)], h[length(h)]]
    empty!(r0)
    push!(r0, V[end], m[end], n[end], h[end])
    println(V[end], "  V[length(V)]  ", sizeof(V), " sizeof(V)")
    println(V[end], " ", m[end]," ", n[end]," ", h[end])
    println(r0, " r0 after push")
    push!(I_appArray, I_app)
    println(I_appArray, " I_appArray")
    println(sizeof(V), " Sizeof(V)")
    global I_app += 0.05
end;
#push!(r0, -58.7085, 0.0953, 0.000913, 0.3662)
using Plots
pyplot()
plot!(I_appArray, AverageFrequency, title = "Bifurcation diagram.",
    color = "blue", xlabel = "I app, nA",
    ylabel = ("ω = 1/T, sec⁠^-1"),
    label = ("Stable Focus"))
gui()
empty!(I_appArray)
empty!(AverageFrequency)
#empty!(r0)


#r0 = [-58.7085; 0.0953; 0.000913; 0.3662];
while I_app >= -0.0005
    t = 0:dt:tf;
    println(r0, " r0")
    (t, r) = ode45(f, r0, t)
    V = map(v -> v[1], r);
    m = map(v -> v[2], r);
    n = map(v -> v[3], r);
    h = map(v -> v[4], r);
    freq = frequency(V, t)

    length(freq) == 0 ? push!(AverageFrequency, 0) : push!(AverageFrequency, sum(freq)/length(freq))
    println(AverageFrequency, " AverageFrequency")
    #r0 = [V[length(V)], m[length(m)], n[length(n)], h[length(h)]]
    empty!(r0)
    push!(r0, V[end], m[end], n[end], h[end])
    println(V[end], "  V[length(V)]  ", length(V), " sizeof(V)")
    println(V[end], " ", m[end]," ", n[end]," ", h[end])
    println(r0, " r0 after push")
    push!(I_appArray, I_app)
    println(I_appArray, " I_appArray")
    println(length(V), " Sizeof(V)")
    global I_app -= 0.05
end;

#using Plots
#using PyPlot
#PyPlot.ioff() #pyplot() # Choose a backend
#PyPlot.plot(t, V, color="red")
#PyPlot.xlabel("t(sec)")
#PyPlot.ylabel("V(t)")
#PyPlot.xlim(0,3) # This will plot to the plot pane
#show()

#using Plots
#pyplot()
plot!(I_appArray, AverageFrequency, title = "Bifurcation diagram.",
    color = "red", xlabel = "I app, nA",
    ylabel = ("ω = 1/T, sec⁠^-1"),
    label = ("Limit Cicle"))
gui()
#plot(n, V)
