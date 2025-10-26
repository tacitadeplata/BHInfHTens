# Paper figure: only Y1 INSTABILITY
# This version corrects a missing 2*mean(pe) in Y2. Results from Y2 are still not the correct ones. Y2 will not be included in the paper. 

#using CSV, DataFrames, Plots, DelimitedFiles, 
using Interpolations
using LaTeXStrings
using Plots
using DelimitedFiles
using DataFrames
using Latexify
using CSV
using Statistics

#Pkg.add("JLD")

function figHLCDM_instability()#( modelname, dim, Y)

PLOTFIGS = 0
#SAVEERRORS = 0

  plot_font = "Helvetica"
  default(fontfamily=plot_font, framestyle=:semi, xtickfont=font(10), ytickfont=font(10), guidefont=font(18), grid=false,xlims=(0,25) )
#        linewidth=2, #label=nothing,
  notaA = string.("a") 
  notaB = string.("b") 
  notaC = string.("c") 
  notaD = string.("d") 
  notaE = string.("e") 
  notaF = string.("f") 
  notaG = string.("g") 
  notaH = string.("h") 


function relativeerror(ph,eph,pe,epe,p,y)
# both p is an arrays
# Posterior error is estimated as the Root mean square deviation (RMSD) between the posterior and the real data

  ep =  sqrt( mean( (sqrt.(p) .- sqrt.(y)).^2) )# this is an scalar
#  println(" ep : ",ep)
#  println( "(eph/(2*ph)) : ", (eph/(2*ph)) )
#  println( " ep ./ p :  ", mean(ep ./ p))
  return  (eph/(2*ph)) .+ (0.5 .* mean(ep ./ p));
end

function bigY1(post, deltaDD, ph, pe)
    return sqrt.(mean(ph)*  mean(eachcol(post))   ) ./ (sqrt.(deltaDD)).^(1/2*mean(pe));
end

function bigY2(post, deltaDD, hp, he)
    return (mean(ph).^(1 /2*mean(pe))) *sqrt.(mean(eachcol(post))) ./ (sqrt.(deltaDD)).^(1/2*mean(pe));
end

#Reading data from the casdade generated with the motif rules
# Loading dimension data after 27 iterations
  filename = "/Users/jlc/EnElaboracion/Cascading_ReservoirNN/prog/Susskind/dims_ReglaCascada_iter27.txt" 
#  filename = "/dims_ReglaCascada_iter27.txt" 
  df = CSV.read(filename, delim=' ', DataFrame,  header=1);

# Identifica cada dimension
  grainsT = df.grainsT;
#coverage = df.coverage
  Dcg = df.Dcg;
  MotifoverG = df.MotifoverG;
  D = df.D;
  n = size(df)[1]

  nR=300000

########################### Data Dcg ############################################
  dd = "Dcg"
  println("Dimension : ",dd)
# Calculo de las diferencias en Dcg
  dDcg = zeros(Float64,n-1)
  for i in 1:n-1
    dDcg[i] = Dcg[i+1] - Dcg[i]
  end
#Cascada directa D/dD = dD/Dinversa
  IdDD = Dcg[1:n-1] ./ dDcg[1:n-1]
#  println("El 1er elemento de IdDD diverge")
#  println("Poping out  1st element")
  popfirst!(IdDD)
  println("Now the IdDD SERIES RUNS FROM x=2")
# generamos x con la misma longuitud que IdDD
  xdDD = collect(2:1:size(IdDD)[1]+1);
  ydDD = IdDD.^2;

  mC="LCDM"
  Cname="chainLCDM"
  mY="Y1"
  println("Esquema : ",mY)
  input1 = "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/LCDM/summarystats-$Cname-$nR$mY.jls"
  input2= "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/LCDM/mean-posterior$mC-$Cname-$nR$mY.csv"
  posteriorM = readdlm(input2, ',', Float64)
  df = CSV.read(input1, delim=',', DataFrame,header=1)
  ph = df.mean[2];
  pe = df.mean[3];
  eph = df.std[2];   epe = df.std[3];  #error_posteriorM = df.mcse[1];

  bigYD = sqrt.(ph.*posteriorM) ./ (sqrt.(ydDD)).^(1 ./2*pe);

# Extrapolo la función inversa 
#Native Interpolations.jl doesn't provide built-in confidence intervals

  ribigYD = reverse(inv.(bigYD))# INVERTED Dcg DATA
  xx=ribigYD[:,1]
  itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
  etp = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation
  h0 = [ 1 ./etp[26]]
  println(etp[26])
  p11 = Plots.plot(26:-1:-26, 1 ./etp[1:53],linestyle=:dash,linewidth=1,seriescolor=:darkblue)
        Plots.scatter!(2:1:26,reverse(1 ./etp[1:25] ) ,size=(300,300),legends=:none,color=1, shape=:+, markersize = 5)
        Plots.plot!(xdDD,bigYD,tick_direction=:out,color=:red,xlims=(-26,27),linewidth=2)      
        Plots.scatter!([1], [h0] )
        Plots.yticks!([10:20:70;], ["10","30","50", "70"])
        Plots.xlabel!(L"\tau")
        Plots.ylabel!(L"\Upsilon(\tau)")
        asymp = string.(round.(1/etp[26]; digits=2))
        annotate!(-28, 1/etp[26]+1, Plots.text(asymp, :red, :right, 14))
        annotate!(25, 15, Plots.text(notaA,:right, 22, "Helvetica Bold")) 
        vline!([1],seriescolor=:grey)      
        hline!([1/etp[26]],seriescolor=:grey)

#  ELCDMDcgY1 = relativeerror(ph,eph,pe,epe, posteriorM, ydDD)*h0

    
########################### Data D motif  ############################################
# Diferencias sucesivas D (motif)
  DiffD = zeros(Float64,n-1)
  nDiff = length(DiffD)
  for i in 1:nDiff
    DiffD[i] = D[i+1] - D[i]
  end
# Storing Inverse of derivative of Dmotif over Dmotif: IdDDmotif
  IdDDmotifs = D[1:nDiff] ./ DiffD[1:nDiff];
  popfirst!(IdDDmotifs)
#  println("Now the IdDDmotifs SERIES RUNS FROM x=2")
#  println( size(IdDDmotifs)[1] )
  xD = collect(2:1:size(IdDDmotifs)[1]+1);
  yD = IdDDmotifs.^2;


  mC="LCDM"
  Cname="chainmotifLCDM"
  mY="Y1"
  println("Esquema : ",mY)
  input1 = "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/LCDM/summarystats-$Cname-$nR$mY.jls"
  input2= "/Users/jlc/EnElaboracion/Population_LeonardSusskind/BHInfTurbHTens/Turing/data/LCDM/mean-posterior$mC-$Cname-$nR$mY.csv"
  posteriorM = readdlm(input2, ',', Float64)
  df = CSV.read(input1, delim=',', DataFrame,header=1)
  ph = df.mean[2]
  pe = df.mean[3]
  eph = df.std[2];   epe = df.std[3];  #error_posteriorM = df.mcse[1];

# Extrapolo la función inversa 
  bigYmotif = sqrt.(ph*posteriorM) ./ (sqrt.(yD)).^(1/2*pe);
  ribigYmotif = reverse(inv.(bigYmotif)) # INVERTED Dmotif DATA
  xx=ribigYmotif[:,1]

  itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
  etp = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation
  println(1.0 ./etp[26])
  h0 = [ 1.0 ./ etp[26]]
  #  p22 = Plots.plot(26:-1:-14, 1 ./etp[1:41],linestyle=:dash,linewidth=1,seriescolor=:darkblue)
  p12 = Plots.plot(26:-1:-26, 1 ./etp[1:53],linestyle=:dash,linewidth=1,seriescolor=:darkblue)      
        Plots.scatter!(2:1:26,reverse(1 ./etp[1:25] ) ,size=(300,300),legends=:none,color=1, xlimits=(-27,27), shape=:+, markersize = 5)
        Plots.plot!(xD,bigYmotif,tick_direction=:out,color=:red,linewidth=2)       
        Plots.scatter!([1], [h0] )
        Plots.yticks!([-200:200:200;], ["-200","0","200"])
        Plots.xlabel!(L"\tau")
        Plots.ylabel!(L"\Upsilon(\tau)")
        asymp = string.(round.(1 ./etp[26]; digits=2))
        annotate!(-28, 1 ./etp[26], Plots.text(asymp, :red, :right, 14))
        annotate!(25, -170, Plots.text(notaB,:right, 22, "Helvetica Bold")) 
        vline!([1],seriescolor=:grey)      
        hline!([1/etp[26]],seriescolor=:grey)     

##################### SIN INVERSA
  rbigYD = reverse(bigYD)
  xx=rbigYD[:,1]
  itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
  etp = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation
  h0 = [etp[26]]
  println(etp[26])
# Already reversing x
  p21 = Plots.plot(26:-1:-26, etp[1:53], xlims=(-27,27), linestyle=:dash,linewidth=1,seriescolor=:darkblue)
        Plots.scatter!(2:1:26,reverse(etp[1:25] ) ,size=(300,300),legends=:none,color=1, shape=:+, markersize = 5)
        Plots.plot!(xdDD,bigYD,tick_direction=:out,color=:red,xlims=(-26,27),linewidth=2)      
        Plots.scatter!([1], [h0] )
        Plots.xlabel!(L"\tau")
        Plots.ylabel!(L"\Upsilon(\tau)")
        asymp = string.(round.(etp[26]; digits=2))
        annotate!(-27, etp[26] + 1, Plots.text(asymp, :red, :right, 14))
        annotate!(25, -260, Plots.text(notaA,:right, 22, "Helvetica Bold")) 
        vline!([1],seriescolor=:grey)      
        hline!([etp[26]],seriescolor=:grey)


    ########################## SIN INVERSA

  ribigYmotif = reverse(bigYmotif)
  xx=ribigYmotif[:,1]

  itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
  etp = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation
  println(etp[26])
  h0 = [etp[26]]
  p22 = Plots.plot(26:-1:-26, etp[1:53],linestyle=:dash,linewidth=1,seriescolor=:darkblue)      
        Plots.scatter!(2:1:26,reverse(etp[1:25] ) ,size=(300,300),legends=:none,color=1, xlimits=(-27,27), shape=:+, markersize = 5)
        Plots.plot!(xD,bigYmotif,tick_direction=:out,color=:red,linewidth=2)       
        Plots.scatter!([1], [h0] )
        Plots.xlabel!(L"\tau")
        Plots.ylabel!(L"\Upsilon(\tau)")
        asymp = string.(round.(etp[26]; digits=2))
        annotate!(-28, etp[26], Plots.text(asymp, :red, :right, 14))
        annotate!(25, 198, Plots.text(notaB,:right, 22, "Helvetica Bold")) 
        vline!([1],seriescolor=:grey)      
        hline!([etp[26]],seriescolor=:grey)     

    ##################################
    ## valores estandarizados: cg
    # Calculo de H0 a normalizando la desviacion de la media con la desviación estandard 
    MeanbigYD = mean(bigYD)
    stdbigYD = std(bigYD)
    snormalbigYD = (bigYD .- MeanbigYD) ./ stdbigYD 

    rsnormalbigYD = reverse(snormalbigYD)
    xx = rsnormalbigYD[:,1]

    itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
    etpsnD = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation  
    println("H_0 : ", (etpsnD[26] * stdbigYD) +  MeanbigYD ) # SIN INVERTIR
    h0nD = (etpsnD[26] * stdbigYD) +  MeanbigYD  # SIN INVERTIR
    println(" h0n : ", h0nD)
    p31 = Plots.plot(26:-1:-26, (etpsnD[1:53].* stdbigYD) .+ MeanbigYD, xlimits=(-27,27),linestyle=:dash,linewidth=1,seriescolor=:darkblue,legends=:none)
    Plots.scatter!(2:26, (stdbigYD .*reverse(etpsnD[1:25])) .+ MeanbigYD, legends=:none, color=1, xlimits=(-27,27), shape=:+, markersize = 7)
   Plots.plot!(xD, bigYD,tick_direction=:out,color=:red,linewidth=2)
   Plots.scatter!([1], [h0nD] )
        Plots.xlabel!(L"\tau")
        Plots.ylabel!(L"\Upsilon(\tau)")
  asymp = string.(round.(h0nD; digits=2))
  annotate!(-28, [h0nD], Plots.text(asymp, :red, :right, 14))
  annotate!(25, -260, Plots.text(notaA,:right, 22, "Helvetica Bold")) 
      vline!([1],seriescolor=:grey)      
      hline!([h0nD],seriescolor=:grey)  
    

    ###############################
## valores estandarizados: motif
    # Calculo de H0 a normalizando la desviacion de la media con la desviación estandard : (bigYmotif .- MeanbigYmotif) ./ stdbigYmotif
    MeanbigYmotif = mean(bigYmotif)
    stdbigYmotif = std(bigYmotif)
    snormalbigYmotif = (bigYmotif .- MeanbigYmotif) ./ stdbigYmotif 

    rsnormalbigYmotif = reverse(snormalbigYmotif)
    xx = rsnormalbigYmotif[:,1]

    itp = Interpolations.interpolate(xx, BSpline(Quadratic(Periodic(OnCell()))   ))
    etpsn = Interpolations.extrapolate( itp, Interpolations.Line() ) #Line - linear extrapolation  
    println("H_0 : ", (etpsn[26] * stdbigYmotif) +  MeanbigYmotif ) # SIN INVERTIR
    h0n = (etpsn[26] * stdbigYmotif) +  MeanbigYmotif  # SIN INVERTIR

    println(" h0n : ", h0n)
#OJO desde el 1er valor: 1:1:25 
    p32 = Plots.plot(26:-1:-26, (etpsn[1:53].* stdbigYmotif) .+ MeanbigYmotif, xlimits=(-27,27),linestyle=:dash,linewidth=1,seriescolor=:darkblue,legends=:none)
    Plots.scatter!(2:26, (stdbigYmotif .*reverse(etpsn[1:25])) .+ MeanbigYmotif,legends=:none, color=1, xlimits=(-27,27), shape=:+, markersize = 7)
   Plots.plot!(xD, bigYmotif,tick_direction=:out,color=:red,linewidth=2)
   Plots.scatter!([1], [h0n] )
        Plots.xlabel!(L"\tau")
        Plots.ylabel!(L"\Upsilon(\tau)")
  asymp = string.(round.(h0n; digits=2))
  annotate!(-28, [h0n], Plots.text(asymp, :red, :right, 14))
  annotate!(25, 198, Plots.text(notaB,:right, 22, "Helvetica Bold")) 
      vline!([1],seriescolor=:grey)      
      hline!([h0n],seriescolor=:grey)  
    
    

    ###########################
    pT = Plots.plot(p11, p12,
                    p21, p22,
                    p31, p32, 
                            layout=(3,2),size=(700,900),# primary =:false, 
        tickfontsize=14,tick_direction=:out,
        top_margin=[0mm 0mm], 
        left_margin=[7mm 6mm],
        right_margin=[0mm 0mm], 
        bottom_margin=[2mm -4mm])#, markersize = 7)

display(pT)


if PLOTFIGS == 1
  savefig("./HLCDM_figY1_instability.png")
end


print("SALIENDO")
end
