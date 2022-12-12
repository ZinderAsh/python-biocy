import sys

def main(filename):
    f = open(filename, "r")

    segments = 0
    while line := f.readline():
        if line[:1] == 'P':
            read_path(line)
            segments += 1
    print(segments)

def read_path(line):
    line = line.split("\t")
    nodes = line[2].split(",")
    cur = int(nodes[0].replace("+", ""))
    variants = 0
    reference = 0
    for i in range(len(nodes)):
        num = int(nodes[i].replace("+", ""))
        if num != cur + 1:
            print("Skipped", num)
            variants += 1
        else:
            reference += 1
        cur = num

    print("variants:", variants)
    print("reference:", reference)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Requires one argument: filename")
    main(sys.argv[1])
