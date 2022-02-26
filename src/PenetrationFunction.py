import sys


eta = 0.1
theta = 1
constant = int(sys.argv[1])
lam = 0.1
ME = sys.argv[2]
ME_effect_size = 2 #1 with Cc and 2 with CC. we add it to what is multiplied by eta
ME_effect = 1 + ME_effect_size*lam
if __name__ == "__main__":
    risk = None
    if ME == "false":
        risk = eta*(1 + constant*theta)
    if ME == "true":
        risk = eta*(1 + constant*theta + ME_effect)
    print(risk)
