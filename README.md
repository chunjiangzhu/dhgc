### Summary

This is the code for Communication Efficient Distributed Hypergraph Clustering, SIGIR 2021.
It has been tested on Ubuntu 16.04.7 LTS and Mac OS Big Sur Version 11.1.


### Enviroment Set-up
#### install Python packages
```
$ pip install numpy, scipy, pandas, h5py, networks, seaborn
```
#### install Julia packages
```
$ julia
julia> using Pkg
julia> Pkg.add("Laplacians")
julia> Pkg.add("LinearAlgebra")
julia> Pkg.add("MAT")
```
Replace the ~/.julia/packages/Laplacians/K6Pgk/src/solverInterface.jl with solverInterface.jl in GitHub  


#### compile c\#
Use the precompiled .exe file in ./HyperReplica/HyperReplica/bin/Release/HyperReplica.exe, or compile by yourself with Visuo Studio. 

#### install mono 
mono is used to call the compiled .exe file for conductance calculation.   
  
For linux:  
As described in https://www.mono-project.com/download/stable/#download-lin  
```
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
sudo apt install apt-transport-https ca-certificates
echo "deb https://download.mono-project.com/repo/ubuntu stable-xenial main" | sudo tee /etc/apt/sources.list.d/mono-official-stable.list
sudo apt update
```

For Mac:    
https://www.mono-project.com/download/stable/#download-mac  

### Run experiment
```
python pagerank.py
```
You can edit the dataset, num_sites and tune parameter c.  

### Acknowledgement
The code for conductance calculation is from [Hypergraph_clustering_based_on_PageRank](https://github.com/atsushi-miyauchi/Hypergraph_clustering_based_on_PageRank)
