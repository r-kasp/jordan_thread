FLAGS = -fsanitize=address -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

start: a.out
%.o: main%.o
	g++ $(FLAGS) $< -o -O3 $@
%.o: %.cpp 
	g++ -c $(FLAGS) $< -o $@

main.o: main.cpp
matrix.o: matrix.cpp matrix.h
thread.o: thread.cpp thread.h

a.out: main.o matrix.o thread.o
	g++ main.o matrix.o thread.o -lasan -o a.out

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.exe
