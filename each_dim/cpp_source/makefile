eigen_dir ?= C:/eigen/

# 通常 ver
pci_cluster_dim_ex: ward.h ward.cpp SI_step.h SI_step.cpp main_ex_step.cpp
	g++ -o ../pci_cluster_dim_ex.out ward.cpp SI_step.cpp main_ex_step.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG
	g++ -o ../pci_cluster_dim_ex_parallel.out ward.cpp SI_step.cpp main_ex_step.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG -fopenmp