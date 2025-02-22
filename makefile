all: res1

res1: m_main.o matrix_operations.o eigen_alg.o
	g++ m_main.o matrix_operations.o eigen_alg.o -o res1

m_main.o: m_main.cpp
	g++ -c m_main.cpp

matrix_operations.o: matrix_operations.cpp
	g++ -c matrix_operations.cpp

eigen_alg.o: eigen_alg.cpp
	g++ -c eigen_alg.cpp