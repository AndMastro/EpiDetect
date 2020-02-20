import sys


eta = 0.1
theta = 2
constant = int(sys.argv[1])
if __name__ == "__main__":
    risk = eta*(1 + constant*theta)
    print(risk)
