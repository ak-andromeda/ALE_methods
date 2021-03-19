import os, re, sys

def extract_logl(rec_file):
    inh = open(rec_file)
    for line in inh:
        if line.startswith(">logl"):
            fields = re.split("\s+", line.rstrip())
            likelihood = fields[1]
    inh.close()
    return likelihood

def check_have_reconciliations(filename, dir_list):

    have_them = 1
    for dir in dir_list:
        if os.path.exists(dir + "/" + filename):
            continue
        else:
            have_them = 0
    return have_them

def make_consel_file(dirs):

    if len(dirs) < 2:
        print("""Usage: write_consel_file.py <space-delimited list of
                 directories containing identically named
                 reconciliation files to compare>""")
        quit()

    score_num = 0
    to_do = [file for file in os.listdir(dirs[0]+"/") if file.endswith("ml_rec")]
    for file in to_do:
        have_them = check_have_reconciliations(file, dirs)
        if have_them == 1:
            for dir in dirs:
                tag = "#" + dir
                lnl = extract_logl(dir + "/" + file)
                if tag in trees:
                    trees[tag] = trees[tag] + str(lnl) + " "
                    summed_lnls[tag] = summed_lnls[tag] + float(lnl)
                else:
                    trees[tag] = str(lnl) + " "
                    summed_lnls[tag] = float(lnl)
            score_num += 1
    scores = ''
    print(str(len(dirs)) + " " + str(score_num))
    for dir in dirs:
        tag = "#" + dir
        print(trees[tag] + "_" + str(summed_lnls[tag]) + "\n" + trees[tag])

if __name__ == "__main__":
    dirs = sys.argv[1:]
    trees = {}
    summed_lnls = {}
    make_consel_file(dirs)
