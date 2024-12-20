import argparse

from classic import fitting_align, local_align
from space_efficient import space_efficient_fitting_align, space_efficient_local_align


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--type", type=str, default="classic", choices=["space_efficient", "classic"])
    parser.add_argument("--mode", type=str, default="local", choices=["fitting", "local"])
    parser.add_argument("--data", type=str, default="8")
    args = parser.parse_args()

    print("[INFO] You're running ** {} {} algorithm ** ".format(args.type, args.mode))
    
    with open(".\\data\\{}.txt".format(args.data), "r") as fp:
        content = fp.readlines()
        v = content[0][:-1]
        w = content[1]
    
    keys = ['A', 'C', 'T', 'G', '-']
    delta = {}
    for i in range(len(keys)):
        delta[keys[i]] = {k : v for (k,v) in zip(keys, [1 if keys[i] == keys[j]  else -1 for j in range(len(keys))])}
    
    print("[INFO] Here we are using N.O.{} data sample. (You can add more) ".format(args.data))
    print("[INFO] Input data: ")
    print("v: {} ".format(v))
    print("w: {} ".format(w))
    print("===================================================\n")

    if args.type == "space_efficient":
        if args.mode == "fitting":
            score, aligned_v, aligned_w = space_efficient_fitting_align(v, w, delta)
        elif args.mode == "local":
            score, aligned_v, aligned_w = space_efficient_local_align(v, w, delta)
        else:
            raise NotImplementedError()
    elif args.type == "classic":
        if args.mode == "fitting":
            score, alignment = fitting_align(v, w, delta)
            aligned_v, aligned_w = alignment.split("\n")
        elif args.mode == "local":
            score, alignment = local_align(v, w, delta)
            aligned_v, aligned_w = alignment.split("\n")
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()

    print("\n===================================================")
    print("[INFO] Results: ")
    print("[INFO] Alignment Score = {} ".format(score))
    print("[INFO] v's alignment: {} ".format(aligned_v))
    print("[INFO] w's alignment: {} ".format(aligned_w))