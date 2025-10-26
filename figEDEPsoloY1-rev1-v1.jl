#Code and data from the article "Black hole information turbulence and the Hubble tension"
#Copyright (C) 2025  Juan Luis Cabrera Fernández

#This project is licensed under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International - Creative Commons licence (CC BY-NC-SA 4.0). 
# https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en


#using CSV, DataFrames, Plots, DelimitedFiles, 
using Interpolations
using LaTeXStrings
using Plots
using DelimitedFiles
using DataFrames
using Latexify
using CSV
using Statistics

function figEDEP()
# Change if you want to plot the figure or save the error file
PLOTFIGS = 1
SAVEERRORS = 0
CAPITALS=0
ITALIC=0


  plot_font = "Helvetica"
  default(fontfamily=plot_font, framestyle=:semi, xtickfont=font(10), ytickfont=font(10), guidefont=font(18), grid=false,xlims=(0,25) ) 


if CAPITALS == 1
    notaA = string.("A") 
    notaB = string.("B") 
    notaC = string.("C") 
    notaD = string.("D") 
    notaE = string.("E") 
    notaF = string.("F") 
    notaG = string.("G") 
    notaH = string.("H") 
elseif ITALIC==1
#    plot_font = "Arial"
    notaA = string.("(a)") 
    notaB = string.("(b)") 
    notaC = string.("(c)") 
    notaD = string.("(d)") 
    notaE = string.("(e)") 
    notaF = string.("(f)") 
    notaG = string.("(g)") 
    notaH = string.("(h)") 
else
    plot_font = "Helvetica"
    notaA = string.("a") 
    notaB = string.("b") 
    notaC = string.("c") 
    notaD = string.("d") 
    notaE = string.("e") 
    notaF = string.("f") 
    notaG = string.("g") 
    notaH = string.("h") 
end

# This is is the error propagation for bigY (both for D and Dmotif).
# The error in posteriorM is estimated as its distance to y:deltaDD/D
# To produce the error bars to each point of Y(tau) it must be multiplied by Y()    
function evolving_error(errY,ph,eph,pe,epe,p,y)
# Parameters are (array_eerror, ph,eph,pe,epe, posteriorM, yD)
        # Parameters are (array_eerror, ph,eph,pe,epe, posteriorM, D)
        println("y must be D with no sqrt (!)")
      errY = abs.( (sqrt.(p) .- y)) # the difference between each pair at fixed time
#  println("  evolving_error(): eerror : ",eerror)
#  println( "(eph/(2*ph)) : ", (eph/(2*ph)) )
#  p es posteriorM . This expresion is obtained propagating the error of Upsilon
  return  (eph/(2*ph)) .+ (0.5 .* errY ./ p);
end
    
#Reading data from the casdade generated with the motif rules
# Loading dimension data after 27 iterations
  filename = "./dims_ReglaCascada_iter27.txt" 
  df = CSV.read(filename, delim=' ', DataFrame,  header=1);

# Identifica cada dimension. Identifies each dimension
  grainsT = df.grainsT;
  Dcg = df.Dcg;
  MotifoverG = df.MotifoverG;
  D = df.D;
  n = size(df)[1]

  nR=300000

########################### Data Dcg ############################################
  dd = "Dcg"
  println("Dimension : ",dd)
# Calculo de las diferencias en Dcg, Dcg differences
  dDcg = zeros(Float64,n-1)
  for i in 1:n-1
    dDcg[i] = Dcg[i+1] - Dcg[i]
  end
#Cascada directa D/dD = dD/Dinversa. Direct cascade D/dD = dD/Dinverse
  IdDD = Dcg[1:n-1] ./ dDcg[1:n-1]
#  println("El 1er elemento de IdDD diverge")
#  println("Poping out  1st element")
  popfirst!(IdDD)
  println("Now the IdDD SERIES RUNS FROM x=2")
# generamos x con la misma longuitud que IdDD. x with same length than IdDD
  xdDD = collect(2:1:size(IdDD)[1]+1);
  ydDD = IdDD.^2;
#Loads summary and posterior info from Turing fitting. 
#Dcg + Y1
  mC="EDEP"
  Cname="chainEDEPD"
  mY="Y1"
  println("Esquema : ",mY)
  input1 = "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/EDEP/summarystats-$Cname-$nR$mY.jls"
  input2= "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/EDEP/mean-posterior$mC-$Cname-$nR$mY.csv"
  posteriorM = readdlm(input2, ',', Float64)
  df = CSV.read(input1, delim=',', DataFrame,header=1)
  ph = df.mean[2]
  pe = df.mean[3]
    eph = df.std[2];
    epe = df.std[3];  #error_posteriorM = df.mcse[1];

    println(" ")
      println("Posterior Parameters and their errors: ")
  println(" ph ", ph, " eph ", eph, " pe ", pe, " epe ", epe )

  p11 = scatter(xdDD, sqrt.(ydDD); label="", color=1, legend=:false,tick_direction=:out,
                         xlimits=(0,27), shape=:+, markersize = 5)
           Plots.plot!(xdDD,sqrt.(posteriorM); label="", color=:red , grid=false, size=(400,400),linewidth=2)   
           Plots.ylabel!(L"\left(\dfrac{\Delta D_{cg}}{D_{cg}}\right)_{I}")    
           Plots.xlabel!(" ")
           Plots.yticks!([0:2200:2200;], ["0", "2200"])
if ITALIC == 1
        annotate!(6, 2200, Plots.text(notaA,:right, 20, "Times Italic")) 
else
        annotate!(5, 2200, Plots.text(notaA,:right, 20, "Helvetica Bold")) 
end

#          display(p11)

          bigYD = sqrt.(ph.*posteriorM) ./ (sqrt.(ydDD)).^(1 ./2*pe);
# Extrapolo la función inversa. Extrapolating the inverse function
#Native Interpolations.jl doesn't provide built-in confidence intervals
          ribigYD = reverse(inv.(bigYD))
          xx=ribigYD[:,1]
          itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
          etp = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation
          h0 =  1 ./etp[26]

  errY = zeros(Float64,n-1)
  errY = evolving_error(errY,ph,eph,pe,epe, posteriorM, IdDD )
  errY = errY .* bigYD
    
# The error in the extrapolation is estimated as the mean of errY in the large fluctuation domain, i.e, errY[1:4]
 println("The error in the extrapolation is estimated as the mean of errY in the large fluctuation domain, i.e, errY[1:4]")
    err_extra = mean(errY[1:4])
#          err_extra = maximum(errY)

    println("Dcd : h_0 = $h0 ± $err_extra")
    p12 = Plots.scatter(2:1:26,reverse(1 ./etp[1:25] ), yerr = errY, markerstrokecolor=1,
                        size=(300,300),legends=:none,color=1, shape=:+, markersize = 5)
     Plots.plot!(xdDD,bigYD,tick_direction=:out,color=:red,xlims=(0,27),linewidth=2)
      Plots.scatter!([1], [h0], yerr= err_extra, markerstrokecolor=:green, markerstrokewidth=1 )
     Plots.ylabel!(L"\Upsilon(z)")
     Plots.xaxis!(minorticks=5)

     asymp = string.(round.(1/etp[26]; digits=2))
     annotate!(-1, 1/etp[26], Plots.text(asymp, :red, :right, 14))
if ITALIC == 1
      annotate!(25, 55, Plots.text(notaB,:right, 20,  "Times Italic"))
else
      annotate!(25, 55, Plots.text(notaB,:right, 20, "Helvetica Bold")) 
end

              Plots.yticks!([50:20:70;], ["50", "70"])
              Plots.xticks!([0:5:25;], ["0","5","10","15","20","25"])
#          display(p12)

#         println("showing pYD")
########################### Data D motif  ############################################
#Dmotif + Y1
# Diferencias sucesivas -  D (motif) -  Differences 
  DiffD = zeros(Float64,n-1)
  nDiff = length(DiffD)
  for i in 1:nDiff
    DiffD[i] = D[i+1] - D[i]
  end
# Storing Inverse of derivative of Dmotif over Dmotif: IdDDmotif
  IdDDmotifs = D[1:nDiff] ./ DiffD[1:nDiff];
  popfirst!(IdDDmotifs)
  println("Now the IdDDmotifs SERIES RUNS FROM x=2")
  println( size(IdDDmotifs)[1] )
  xD = collect(2:1:size(IdDDmotifs)[1]+1);
  yD = IdDDmotifs.^2;


  mC="EDEP"
  Cname="chainDheEDEP"
  mY="Y1"
  println("Esquema : ",mY)
  input1 = "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/EDEP/summarystats-$Cname-$nR$mY.jls"
  input2= "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/EDEP/mean-posterior$mC-$Cname-$nR$mY.csv"
  posteriorM = readdlm(input2, ',', Float64)
  df = CSV.read(input1, delim=',', DataFrame,header=1)
  ph = df.mean[2]
  pe = df.mean[3]
    eph = df.std[2];
    epe = df.std[3];  #error_posteriorM = df.mcse[1];

    println(" ")
        println("Posterior Parameters and their errors: ")
  println(" ph ", ph, " eph ", eph, " pe ", pe, " epe ", epe )
# Posterior error is estimated as the sqrt of the square of the distance between the posterior and the real data

  p21 = Plots.scatter(xD, sqrt.(yD ); label="", color=1, legend=:false,tick_direction=:out, shape=:+, markersize = 5)
        Plots.plot!(xD,sqrt.(posteriorM); label="", color=:red , grid=false,size=(400,400), xlimits=(0,27),linewidth=2)      
            Plots.ylabel!(L"\left(\dfrac{\Delta D_m}{D_m}\right)_{I}")
            Plots.xlabel!(L"\tau")
            Plots.yticks!([200:600:800;], ["200", "800"])#,title="pHLCDM3")
if ITALIC == 1
        annotate!(6, 850, Plots.text(notaC,:right, 20, "Times Italic"))
else
        annotate!(5, 850, Plots.text(notaC,:right, 20, "Helvetica Bold")) 
end
# display(p21)

# Extrapolo la función inversa. Extrapolating the inverse function
    bigYmotif = sqrt.(ph*posteriorM) ./ (sqrt.(yD)).^(1/2*pe);
    ribigYmotif = reverse(inv.(bigYmotif))
    xx=ribigYmotif[:,1]

    itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
    etp = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation
    println(1 ./etp[26])
    h0 =  1 ./etp[26]
    
    errY = evolving_error(errY,ph,eph,pe,epe, posteriorM, IdDDmotifs )
    errY = errY .* bigYmotif
    
# The error in the extrapolation is estimated as the mean of errY in the large fluctuation domain, i.e, errY[1:4]
#    err_extra = mean(errY) NO
    err_extra = mean(errY[1:4])
#    println(errY)
    println("Dmotif : h_0 = $h0 ± $err_extra")
    p22 = Plots.scatter(2:1:26,reverse(1 ./etp[1:25] ), yerr=errY, markerstrokecolor=1,
                        size=(300,300),legends=:none,color=1, xlimits=(0,27), shape=:+, markersize = 5)
          Plots.plot!(xD,bigYmotif,tick_direction=:out,color=:red,linewidth=2)       
          Plots.scatter!([1], [h0], yerr=err_extra, markerstrokecolor=:green)#, markerstrokewidth=1,markersize=3 )
          Plots.ylabel!(L"\Upsilon(z)")#,guide_position=:right,guidefontrotation=180)
          Plots.xlabel!(L"\tau")
          xaxis!(minorticks=5)
          asymp = string.(round.(1 ./etp[26]; digits=2))
          annotate!(-1, 1 ./etp[26], Plots.text(asymp, :red, :right, 14))
if ITALIC == 1
  annotate!(25, 36.5, Plots.text(notaD,:right, 20, "Times Italic"))
else
  annotate!(25, 36.5, Plots.text(notaD,:right, 20, "Helvetica Bold")) 
end

Plots.yticks!([25:10:35;], ["25", "35"])

# display(p22)
  println("El elemento 26 a tau=1 : ", h0)
##############################################################################
  pT = Plots.plot(p11, p12, p21, p22,       
        layout=(2,2),size=(700,700), 
        tickfontsize=14,tick_direction=:out,
        top_margin=[0mm 0mm], 
        left_margin=[7mm 6mm],
        right_margin=[0mm 0mm], 
                  bottom_margin=[2mm -4mm])
display(pT)
if PLOTFIGS == 1
#    savefig("./PRLfig10Y1.png")
    savefig("./figEDEPsoloY1-rev1-v1.png")
end

    print("GOING OUT ...")
end
