gcc -Wall -m64 -g -o L1_MIP ./L1_MIP.c  -I/opt/gurobi550/linux64/include/ -L/opt/gurobi550/linux64/lib/ -lgurobi55 -lepanet -lpthread -lm && ./L1_MIP
