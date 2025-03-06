using IntervalArithmetic
using Plots 
using Dates
using StatsBase 
using Base.Threads

function phi_n(phi,T,x, s, n)
    "Calculate the image of phi^{(n)}"
    if n == 0
        return s
    else
        return phi.image(x, phi_n(phi,T,T.image(x), s, n-1))
    end
end

function phi_nBoundK(phi,T,x, s, n,K)
    "Calculate the image of phi^{(n)} bound by K"
    if n == 0
        return s
    else
        return phi.image(x,min(phi_nBoundK(phi,T,T.image(x), s, n-1,K),K))
    end
end

function phi_n_per(phi,orb,len,s,n,i,err)
    "Calculate the image of phi^{(n)} for a periodic orbit"
    if n == 0
        return interval(s,s)
    else
        return phi.image(interval(orb[i+1]-err,orb[i+1]+err), phi_n_per(phi,orb,len, s, n-1,(i+1)%len,err))
    end
end

function compose(f,x,n)
    "Calculate the image of the n-th compound of f"
    if n == 0
        return x
    else
        return compose(f, f(x), n-1)
    end
end

function split_interval(I::Interval)
    mid = (inf(I) + sup(I)) / 2
    return interval(inf(I), mid), interval(mid, sup(I))
end

## Upper bound of LambdaF

function FindK(phi,T,N)
    "find a good constant upper than q"
    K=1/2
    A=[interval((i-1)/2^N,i/2^N) for i in 1:2^N]
    while K<0.9999
        r=Inf
        for NK in 1:12
            ph=maximum([sup(phi_n(phi,T,x, K, NK)) for x in A])
            if K>ph
                r=min(r,ph)
            end
        end
        if r!=Inf
            return r
        end
        K=(1+K)/2
    end
    error("K not find")
end

function upboundlambdaF(phi,T,n,N,m,Q,I)
    "find a good upper bound of lambdaF"
    Indice=[i for i in 1:2^N]
    Nmaj=N
    bound=-Inf
    K=FindK(phi,T,12)
    for z in 1:I
        len=length(Indice)
        l=[interval((i-1)/2^Nmaj,i/2^Nmaj) for i in Indice]
        Tl=[]
        r=[interval(0,0) for i in 1:len]
        infi=[0 for i in 1:len]
        for i in 1:m
            Tl=[T.image(x) for x in l]
            majq=[phi_nBoundK(phi,T,x,interval(K,K),n,K) for x in Tl]
            r+=[phi.logds(l[j],majq[j]) for j in 1:len]
            infi=[min(sup(r[j])/i,infi[j]) for j in 1:len]
            l=Tl
        end
        boundmaj=quantile(infi,1-1/2^Q)
        if boundmaj==0
            return 0
        end
        if bound>boundmaj 
            return max(maximum(infi),bound)
        end
        if z==I
            return maximum(infi)
        end
        bound=boundmaj
        Nmaj=Nmaj+Q
        Indicemaj=Indice[findall(x -> x > bound, infi)].*2^Q
        Indice=reduce(vcat,[Indicemaj.-i for i in 0:2^Q-1])
    end
end

function upboundlambdaF2(phi,T,n,N,m,I,NK)
    "find a good upper bound of lambdaF"
    linit=[interval((i-1)/2^N,i/2^N) for i in 1:2^N]
    K=FindK(phi,T,NK)
    for z in 1:I
        len=length(linit)
        l=linit
        Tl=[]
        r=[interval(0,0) for i in 1:len]
        infi=[Inf for i in 1:len]
        for i in 1:m
            Tl=[T.image(x) for x in l]
            majq=[phi_nBoundK(phi,T,x,interval(K,K),n,K) for x in Tl]
            r+=[phi.logds(l[j],majq[j]) for j in 1:len]
            infi=[min(sup(r[j])/i,infi[j]) for j in 1:len]
            l=Tl
        end
        bound=quantile(infi,0.95)
        if z==I
            return maximum(infi)
        end
        lbis=[]
        for i in 1:len
            if infi[i]>bound
                a,b=split_interval(linit[i])
                push!(lbis,a)
                push!(lbis,b)
            else
                push!(lbis,linit[i])
            end
        end
        linit=lbis
    end
end

function upboundlambdaF3(phi,T,N,I,NK)
    "find a good upper bound of lambdaF"
    linit=[inte(interval((i-1)/2^N,i/2^N),N,1,1000.0) for i in 1:2^N]
    len=length(linit)
    K=FindK(phi,T,NK)
    for z in 1:I
        for x in linit
            if x.new==1
                r=interval(0,0)
                m=div(x.len,2)+3
                n=m
                inter=x.in
                for j in 1:m
                    Tinter=T.image(inter)
                    majq=phi_nBoundK(phi,T,Tinter,interval(K,K),n,K)
                    r+=phi.logds(inter,majq)
                    x.bound=min(sup(r)/j,x.bound)
                    inter=Tinter
                end
            end
        end
        if z==I
            return maximum([x.bound for x in linit])
        end
        bound=quantile([x.bound for x in linit],0.95)
        "lbis=[]
        for x in linit
            if x.bound>bound
                a,b=split_interval(x.in)
                x1=inte(a,x.len+1,1,x.bound)
                x2=inte(b,x.len+1,1,x.bound)
                push!(lbis,x1)
                push!(lbis,x2)
            else
                x=inte(x.in, x.len, 0, x.bound)
                push!(lbis,x)
            end
        end"
        len+=count(x->x.bound>bound,linit)
        lbis=Vector{Any}(undef, len)
        j=1
        for x in linit
            if x.bound>bound
                a,b=split_interval(x.in)
                x.in=a
                x.len+=1
                x.new=1
                lbis[j]=x
                j+=1
                x1=deepcopy(x)
                x1.in=b
                lbis[j]=x1
            else
                x.new=0
                lbis[j]=x
            end
            j+=1
        end
        linit=lbis
    end
end

## Orbits

"finds a representative of each sequence of size n with integer values 
between 0 and d up to rotation and not containing a repeating pattern"

function is_rotation_of_minimal(seq, n)
    for i in 1:n-1
        rotated_seq = vcat(seq[i+1:end], seq[1:i])
        if rotated_seq < seq
            return false
        end
    end
    return true
end

function has_repeating_pattern(seq, n)
    for len in 1:div(n, 2)
        if n % len == 0
            pattern = seq[1:len]
            is_repeating = true
            for i in 1:div(n, len)
                if seq[(i-1)*len+1:i*len] != pattern
                    is_repeating = false
                    break
                end
            end
            if is_repeating
                return true
            end
        end
    end
    return false
end

function generate_minimal_non_repeating_sequences(d::Int, n::Int)
    minimal_sequences = Vector{Vector{Int}}()  
    seq = fill(0, n)  
    while seq != nothing
        if is_rotation_of_minimal(seq, n) && !has_repeating_pattern(seq, n)
            push!(minimal_sequences, copy(seq))  # Ajouter une copie de la séquence à la liste
        end
        seq = next_sequence(seq, d)
    end
    
    return minimal_sequences
end

function next_sequence(seq::Vector{Int}, d::Int)
    n = length(seq)
    for i in n:-1:1
        if seq[i] < d - 1
            seq[i] += 1
            for j in i+1:n
                seq[j] = 0
            end
            return seq
        end
    end
    return nothing
end

function dichotomy(g,z,xm,XM,epsi)
    "find an approximation of g(x)=z with an error less than epsi"
    t=(xm+XM)/2
    while XM-xm>epsi
        if g(t)>z
            XM=t
        else
            xm=t
        end
        t=(xm+XM)/2
    end
    return t
end

function suite_dadique(suite::Vector{Int}, d::Int,i)
    value = [sum(suite[(j+i-k)%i+1] * d^(j-1) for j in 1:i) for k in 1:i]
    return value
end

function logdphiq(phi,orbx,i,N,n,epsi)
    lf=0
    qTx=phi_n_per(phi,orbx,i,0,n,0,epsi)
    for j in i:-1:1
        lf+=inf(phi.logds(interval(orbx[j]-epsi,orbx[j]+epsi),qTx))
        qTx=phi.image(interval(orbx[j]-epsi,orbx[j]+epsi),qTx)
    end
    return lf/i
end

## Upper bound on Lambda_u

function upperboundLu(T,N,m,Q,I)
    "find a good upper bound of lambda_u"
    Indice=[i for i in 1:2^N]
    Nmaj=N
    bound=-Inf
    for z in 1:I
        len=length(Indice)
        if len==0
            return bound
        end
        l=[interval((i-1)/2^Nmaj,i/2^Nmaj) for i in Indice]
        r=[interval(0,0) for i in 1:len]
        infi=[Inf for i in 1:len]
        for i in 1:m
            r+=[log(T.d(x)) for x in l]
            infi=[min(sup(r[j])/i,infi[j]) for j in 1:len]
            l=[T.image(x) for x in l]
        end
        boundmaj=quantile(infi,1-1/2^Q)
        if bound>boundmaj 
            error(bound)
        end
        bound=boundmaj
        if z==I
            return maximum(infi)
        end
        Nmaj=Nmaj+Q
        Indicemaj=Indice[findall(x -> x > bound, infi)].*2^Q
        Indice=reduce(vcat,[Indicemaj.-i for i in 0:2^Q-1])
    end
end

function logdTseq(T,orbx,epsi)
    r=interval(0,0)
    i=length(orbx)
    for j in 1:i
        r+=log(T.d(interval(orbx[j]-epsi,orbx[j]+epsi)))
    end
    return inf(r)/i
end

## Plot Holder regularity

function plotBoundHolderreg(phi,T,premierl,pas,nbl,n,N,m,M,nu,Nu,plt)
    "step 1: find approximation of the periodique orbit"
    d=T.degree
    epsi=1/(2^N)
    Orbite=Vector{Vector{Float64}}()
    push!(Orbite,[0])
    for i in 1:M
        seq=[]
        g=x->compose(T.image,x,i)-x
        if i>1
            seq=generate_minimal_non_repeating_sequences(T.degree,i)
        end
        if i==1 && d>2
            seq=[[j] for j in 1:d-2]
        end
        for x in seq 
            image_orbx=suite_dadique(x,T.degree,i)
            orbx=[]
            xm=0
            xM=1
            for f in image_orbx
                a=dichotomy(g,f,xm,xM,epsi)
                push!(orbx,a)
                Ta=T.image(a)
                err=epsi*T.maxd
                xm=Ta-err
                xM=Ta+err
            end
        push!(Orbite,orbx) 
        end
    end
    "step 2: find Lu"
    lowLu=maximum([logdTseq(T,orbx,epsi) for orbx in Orbite])
    upLu=upperboundLu(T,10,14,1,20)
    "step 3: find Lf"
    upLf=zeros(Float64,nbl)
    upLf2=zeros(Float64,nbl)
    @time begin
    for p in 1:nbl
        l=premierl+(p-1)*pas
        upLf[p]=upboundlambdaF3(phi(l),T,4,160,14)
    end
    end
    @time begin
    for p in 1:nbl
        l=premierl+(p-1)*pas
        upLf2[p]=upboundlambdaF2(phi(l),T,16,4,16,90,14)
    end
    end
    lowLf=zeros(Float64,nbl)
    for p in 1:nbl
        l=premierl+(p-1)*pas
        lowLf[p]=maximum([logdphiq(phi(l),orbx,length(orbx),N,200,epsi) for orbx in Orbite])
    end
    lowReg=upLf./(-upLu)
    lowReg2=upLf2./(-upLu)
    upReg=lowLf./(-lowLu)
    if plt==true
        L=[premierl + i*pas for i in 0:nbl-1]
        p=scatter(L,upReg,color=:blue)
        scatter!(p,L,lowReg2,color=:green)
        scatter!(p,L,lowReg,color=:red)
        times=Dates.format(now(),"yyyy-mm-dd_HH-MM-SS")
        ###savefig(homedir() * "/Bureau/plot/" * "graphe_$times.png")
        display(p)
    end
    return upReg,lowReg2,lowReg
end
    
## Examples of transformations

struct T
    image
    d
    maxd
    mind
    degree
end

T1(e)=T(x->2*x+e*sin(2*pi*x),x->2+2*pi*e*cos(2*pi*x),2+2*pi*e,2-2*pi*e,2)

T2(n)=T(x->n*x,x->n,n,n,n)

## Examples of laws of reproduction

struct phi
    image
    logds
    exp
    maxexp
    maxdexp
    maxdx
end

phi1(l)=phi((x,s)->exp((s - 1) * exp(l + cos(2 * x * π))),(x,s)->log(exp(l + cos(2 * x * π)) * phi1(l).image(x, s)),x->exp(l+cos(2 * x * π)),exp(l+1),2*π*exp(l+1),2*π/exp(1))

phi2(l)=phi((x,s)->exp((s - 1) * exp(l - cos(2 * x * π))),(x,s)->log(exp(l - cos(2 * x * π)) * phi2(l).image(x, s)),x->exp(l-cos(2 * x * π)),exp(l+1),2*π*exp(l+1),2*π/exp(1))


mutable struct inte
    in
    len
    new
    bound
end
## General settings 

N = 14 # 2^N est le nb d'intervalles
n = 14 # le nb d'itérations de notre fonction majorante et minorante
m = 8
M = 8 # taille max des orbites périodiques

Nu=14
nu=8
Mu=8

e=0.1
"e doit etre inférieur à 1/(pi*2) soit environ 0.159"
l=0.6

## Bounds on q 

function upboundq(phi,T,n,N,plt)  # N>n
    "find a good upper bound of q"
    K=FindK(phi,T,N)
    l1=[[interval(i/2^N,(i+1)/2^N),interval(K,K)]  for i in 0:(2^N-1)]
    l2=[phi_n(phi,T,x[1],x[2],n) for x in l1]
    if plt==true
         x=[i/2^N for i in 1:2^N]
        majqm=[sup(i) for i in l2]
        p=plot(x,majqm,title="upboundq",ylim=(-0.01,1.01))
        display(p)
    end
    return l2
end

function lowboundq(phi,T,n,N,plt)  # N>n
    "find a good lower bound of q"
    l1=[[interval(i/2^N,(i+1)/2^N),interval(0,0)]  for i in 0:(2^N-1)]
    l2=[phi_n(phi,T,x[1],x[2],n) for x in l1]
    if plt==true
        x=[i/2^N for i in 1:2^N]
        minqm=[inf(i) for i in l2]
        p=plot(x,minqm,title="lowboundq",ylim=(-0.01,1.01))
        display(p)
    end
    return l2
end

function boundq(phi,T,n,N,plt)  # N>n
    "find good bounds of q"
    K=FindK(phi,T,N)
    l1=[[interval(i/2^N,(i+1)/2^N),interval(K,K)]  for i in 0:(2^N-1)]
    l2=[phi_n(phi,T,x[1],x[2],n) for x in l1]
    l3=[[interval(i/2^N,(i+1)/2^N),interval(0,0)]  for i in 0:(2^N-1)]
    l4=[phi_n(phi,T,x[1],x[2],n) for x in l3]
    if plt==true
        x=[i/2^N for i in 1:2^N]
        majqm=[sup(i) for i in l2]
        minqm=[inf(i) for i in l4]
        plot(x,majqm,label="upboundq",ylim=(-0.01,1.01))
        display(plot!(x,minqm,label="lowboundq"))
    end
    return l2,l4
end

## Lower bound on LambdaF

function lowboundlambdaF(phi,T,n,N,M)
    "find a good lower bound of lambdaF"
    d=T.degree
    lf=inf(phi.logds(0,phi_n(phi,T,interval(0,0),interval(0,0),n)))
    epsi=min(1/(2^N),1/(2*(1+T.maxd)*(T.maxd^M-1)))/10 
    for i in 1:M
        g=x->compose(T.image,x,i)-x
        if i>1
            seq=generate_minimal_non_repeating_sequences(T.degree,i)
        end
        if i==1
            seq=[[j] for j in 1:d-2]
        end
        for x in seq 
            image_orbx=suite_dadique(x,T.degree,i)
            orbx=[]
            xm=0
            xM=1
            for f in image_orbx
                a=dichotomy(g,f,xm,xM,epsi)
                push!(orbx,a)
                Ta=T.image(a)
                err=epsi*T.maxd
                xm=Ta-err
                xM=Ta+err
            end
            lf=max(lf,logdphiq(phi,orbx,i,N,n,epsi))
        end
    end
    return lf
end

