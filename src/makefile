CC=g++ -O3
CFLAGS=-W -Wall -pedantic -g
LDFLAGS=-lgmp

all: test_add_gmp test_add_idd test_mul_gmp test_mul_idd test_pow_gmp test_pow_idd example

test_add_gmp: test_add_gmp.o iddlib.o
	$(CC) -o test_add_gmp monitor.o iddlib.o test_add_gmp.o $(LDFLAGS)
test_add_gmp.o:	test_add_gmp.cpp
	$(CC) -o test_add_gmp.o -c test_add_gmp.cpp $(CDFLAGS)

test_add_idd: test_add_idd.o iddlib.o
	$(CC) -o test_add_idd monitor.o iddlib.o test_add_idd.o $(LDFLAGS)
test_add_idd.o:	test_add_idd.cpp
	$(CC) -o test_add_idd.o -c test_add_idd.cpp $(CDFLAGS)

test_mul_gmp: test_mul_gmp.o iddlib.o
	$(CC) -o test_mul_gmp monitor.o iddlib.o test_mul_gmp.o $(LDFLAGS)
test_mul_gmp.o:	test_mul_gmp.cpp
	$(CC) -o test_mul_gmp.o -c test_mul_gmp.cpp $(CDFLAGS)

test_mul_idd: test_mul_idd.o iddlib.o
	$(CC) -o test_mul_idd monitor.o iddlib.o test_mul_idd.o $(LDFLAGS)
test_mul_idd.o:	test_mul_idd.cpp
	$(CC) -o test_mul_idd.o -c test_mul_idd.cpp $(CDFLAGS)

test_pow_gmp: test_pow_gmp.o iddlib.o
	$(CC) -o test_pow_gmp monitor.o iddlib.o test_pow_gmp.o $(LDFLAGS)
test_pow_gmp.o:	test_pow_gmp.cpp
	$(CC) -o test_pow_gmp.o -c test_pow_gmp.cpp $(CDFLAGS)
	
test_pow_idd: test_pow_idd.o iddlib.o
	$(CC) -o test_pow_idd monitor.o iddlib.o test_pow_idd.o $(LDFLAGS)
test_pow_idd.o:	test_pow_idd.cpp
	$(CC) -o test_pow_idd.o -c test_pow_idd.cpp $(CDFLAGS)


iddlib.o: iddlib.cpp iddlib.h monitor.o
	$(CC) -o iddlib.o -c iddlib.cpp $(CFLAGS)

monitor.o: monitor.cpp monitor.h 
	$(CC) -o monitor.o -c monitor.cpp $(CFLAGS)
	
example.o: example.cpp
	$(CC) -o example.o -c example.cpp $(CFLAGS)
example: example.o monitor.o iddlib.o
	$(CC) -o example monitor.o iddlib.o example.o $(LDFLAGS)

clean:
	rm -f *.c~ *.h~ *.cpp~ *.ini~ *.o makefile~ test_add_gmp test_add_idd test_mul_gmp test_mul_idd test_pow_gmp test_pow_idd example
