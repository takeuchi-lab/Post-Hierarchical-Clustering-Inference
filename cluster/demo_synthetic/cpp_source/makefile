eigen_dir ?= C:/eigen/

# 通常 ver
pci_cluster_synthetic: ward.h ward.cpp SI_step.h SI_step.cpp main_synthetic.cpp
	g++ -o ../pci_cluster_synthetic ward.cpp SI_step.cpp main_synthetic.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG
	g++ -o ../pci_cluster_synthetic_parallel ward.cpp SI_step.cpp main_synthetic.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG -fopenmp