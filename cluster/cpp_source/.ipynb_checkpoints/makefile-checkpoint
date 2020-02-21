eigen_dir ?= C:/eigen/

# 通常 ver
pci_cluster_ex: ward.h ward.cpp SI_step.h SI_step.cpp main_ex.cpp
	g++ -o ../pci_cluster_ex ward.cpp SI_step.cpp main_ex.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG
	g++ -o ../pci_cluster_ex_parallel ward.cpp SI_step.cpp main_ex.cpp -I $(eigen_dir) -DEIGEN_NO_DEBUG -fopenmp