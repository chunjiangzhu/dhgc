### About

This is the code for Communication Efficient Distributed Hypergraph Clustering, SIGIR 2021.
It has been tested on Ubuntu 16.04.7 LTS and Mac OS Big Sur Version 11.1.


### Enviroment Set-up
#### Install Python packages
```
$ pip install numpy, scipy, pandas, h5py, networkx
```
#### Install Julia packages
```
$ julia
julia> using Pkg
julia> Pkg.add("Laplacians")
julia> Pkg.add("LinearAlgebra")
julia> Pkg.add("MAT")
```
Replace the ~/.julia/packages/Laplacians/K6Pgk/src/solverInterface.jl with Laplacians/solverInterface.jl in GitHub  


#### Compile c\#
Use the precompiled .exe file in ./HyperReplica/HyperReplica/bin/Release/HyperReplica.exe, or compile by yourself with Visual Studio. This is for conductance calculation. 

#### Install mono 
mono is used to call the compiled .exe file for conductance calculation.   
  
For linux:  
As described in https://www.mono-project.com/download/stable/#download-lin  
```
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
sudo apt install apt-transport-https ca-certificates
echo "deb https://download.mono-project.com/repo/ubuntu stable-xenial main" | sudo tee /etc/apt/sources.list.d/mono-official-stable.list
sudo apt update
sudo apt install mono-devel
sudo apt install mono-runtime
```

For Mac:    
https://www.mono-project.com/download/stable/#download-mac  

### Datasets
Datasets are from 2 resources:   
[Cornnel](https://www.cs.cornell.edu/~arb/data/) and [Hypergraph Clustering Based on PageRank](https://github.com/atsushi-miyauchi/Hypergraph_clustering_based_on_PageRank) in [KDD'21](https://dl.acm.org/doi/10.1145/3394486.3403248)

### Run experiment
```
python pagerank.py --dataset highschool --num_sites 3 --c 2 
```
You can edit the dataset, num_sites and tune parameter c.  

#### Suggestions
The program includes some disk I/O. If the program was down in the middle of the process, some of the generated files may not be complete. Please remove the ./data/DATASET/tmp folder and then re-run the program.  

### Acknowledgement
The code for conductance calculation is from [Hypergraph_clustering_based_on_PageRank](https://github.com/atsushi-miyauchi/Hypergraph_clustering_based_on_PageRank)
