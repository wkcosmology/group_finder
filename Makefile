run_group_finder: run_group_finder.cpp
	$(CC) -O3 -std=c++17 $^ -o run_group_finder.out -lhippcntl -lhdf5 -lhippio -lstdc++fs -lm
nbr_finder.o: nbr_finder.cpp
	$(CC) -O3 -std=c++17 -c $<
clean:
	rm *.out *.o
