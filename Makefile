all:
	@gcc gm.c main.c -o main
clean:
	@rm -rf *.o *.out *.exe *.bak