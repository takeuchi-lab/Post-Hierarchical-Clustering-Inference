eigen_dir ?= C:/eigen/

# 通常 ver
pci_cluster_FPR: ward.h ward.cpp SI_step.h SI_step.cpp main_FPR.cpp
	g++ -o ../pci_cluster_FPR ward.cpp SI_step.cpp main_FPR.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG
	g++ -o ../pci_cluster_FPR_parallel ward.cpp SI_step.cpp main_FPR.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG -fopenmp