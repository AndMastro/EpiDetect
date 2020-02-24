import sys


eta = 0.1
theta = 0.5
constant = int(sys.argv[1])
if __name__ == "__main__":
    risk = eta*(1 + constant*theta)
    print(risk)
