from unittest import result
from matplotlib.pyplot import text
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
import sys, os
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import time
from ete3 import Tree, TreeStyle, NodeStyle, AttrFace, faces, TextFace
from collections import defaultdict
import shlex, subprocess
from subprocess import CalledProcessError
# from urllib.request import urlopen
# from requests_html import HTMLSession
# from bs4 import BeautifulSoup
# executable_path='C:/Users/ming/Downloads/chromedriver_win32/chromedriver.exe'



def main(args):

    opts = parse_cmdline_params(args[1:])
    PATH = 'https://www.ncbi.nlm.nih.gov/pathogens/'

    options = webdriver.ChromeOptions()
    prefs = {"download.default_directory" : opts.folders}
    options.add_experimental_option("prefs",prefs)
    
    driver = webdriver.Chrome(chrome_options=options)
    # driver.set_window_position(-10000,0)
    driver.get(PATH)
    # find the searchbox
    searchbox = driver.find_element(By.XPATH,'//*[@id="main-isolates-search-field"]')
    # send SRR to search box
    searchbox.send_keys(opts.SRR)
  
    search = waitUntil(driver, By.XPATH,'//*[@id="main-isolates-search"]/button/span')
    search.click()
    # driver.find_element(By.XPATH,'//*[@id="main-isolates-search"]/button/span').click()

    searchTable = waitUntil(driver, By.XPATH,'//*[@id="gridview-1022-record-65"]/tbody/tr/td[11]/div/a')
    # driver.find_element(By.XPATH,'//*[@id="gridview-1022-record-65"]/tbody/tr/td[6]/div/a')
    searchTable.click()

    driver.switch_to.window(driver.window_handles[1])
    print(driver.current_url)


    targetPD = waitUntil(driver, By.XPATH,'//*[@id="gridview-1022-record-78"]/tbody/tr/td[7]/div')
    targetPD = targetPD.text

    download = waitUntil(driver, By.XPATH,'//*[@id="button-1133"]')
    download.click()
    
    newick = waitUntil(driver, By.XPATH,'//*[@id="menuitem-1135-textEl"]')
    newick.click()
    # driver.find_element(By.XPATH,'//*[@id="menuitem-1135-textEl"]').click()
    driver.close()

    treeRead(opts,targetPD)
    


def treeRead(opts, targetPD):
    os.chdir(opts.folders)  
    print(os.getcwd())
    metafile = opts.file
    
    
    t = Tree("export.newick", quoted_node_names=True, format=1)
    
    print(t)
    pdname = find_target(targetPD, t)
    print()
    l = find_close(pdname,t)
    t.prune(l[0:11])
    a = [node.name.split(", ")[-1] for node in t.traverse("postorder") if node.name.split(", ")[-1]]
    l_matching = findSRR(a, metafile)

    #  add SRR number behind each isolate   
    for node in t.traverse("postorder"):
        # Do some analysis on node
        if node.name.split(", ")[-1] in l_matching and l_matching[node.name.split(", ")[-1]] != "NULL":
            node.name += ", "+l_matching[node.name.split(", ")[-1]] 
            # if targetPD in node.name  :
            #     node.add_face(TextFace(l_matching[node.name.split(", ")[-1]] , fsize=8, ftype = "Arial", fgcolor="green"),position="aligned", column=0)
            # else: 
            #     node.add_face(TextFace(l_matching[node.name.split(", ")[-1]] , fsize=8, ftype = "Arial", fgcolor="black"),position="aligned", column=0)
    print(t)

    # t.write(outfile="new_tree.nwk",format= 1)

    style = TreeStyle()
    style.mode = "r" # draw tree in circular mode
    style.show_leaf_name = True
    style.show_branch_length = True
    style.branch_vertical_margin = 5
    style.scale =  30
    style.margin_top = 15
    nst1 = NodeStyle()
    nst1["bgcolor"] = "yellow"
    nst2 = NodeStyle()
    nst2["fgcolor"] = "red"
    nst2["size"] = 10

    for n in t.traverse():
        if n.is_leaf():
        

            if targetPD in n.name  :
                # name_face = AttrFace("name", fsize=12, fgcolor="#009000")
                 # Adds the name face to the image at the preferred position
                # n.add_face_to_node(name_face, column=0)
                # faces.add_face_to_node(nameFace, n, column=0)
                n.set_style(nst1)

            if "clinical" in n.name and targetPD not in n.name:
                n.set_style(nst2)
        
    # add title
    style.title.add_face(TextFace("Close isolates of " + targetPD, fsize=20), column=0)
    t.render("mytree.png", tree_style=style, dpi = 200)

def find_target(name, tree):
    for node in tree.traverse("postorder"):
        # Do some analysis on node
        if name in node.name:
            print(node.name)
            return node.name
        

def find_close(name0,tree):
    
    

    dist_dict = defaultdict(list)
    l = []
    target = tree.search_nodes(name=name0)[0]

    for node in tree.iter_descendants("postorder"):
        # Do some analysis on node
        
        dist = target.get_distance(node)
        # if dist not in dist_dict:
        #     dist_dict[dist] = node.name

        
        dist_dict[dist].append(node.name)
        


    dist_dict=dict(sorted(dist_dict.items()))
    for key in dist_dict:
        if '' not in dist_dict[key]:
            l.extend(dist_dict[key] )
    
    return l

def findSRR(l_PDT, metafile):

    cmd1 = "awk -F '\t' "
    cmd2 = " { print $9, $42 }' " + metafile
    cmd3 = "'" + "$42 == "+ '"' + l_PDT[0] + '"'
    for pdt in l_PDT[1:11]:
        cmd3 += " || "+"$42 == "+ '"' + pdt + '"'

    f = open("matching.txt", "w+")
    args =  cmd1 + cmd3 + cmd2
    final_args= shlex.split(args)
    try:
        proc = subprocess.call(final_args,stdout=f)
        # proc.wait()
        # (stdout, stderr) = proc.communicate()
    except CalledProcessError as err:
        print("Error ocurred: " + err.stderr)

    match = open('matching.txt', 'r')
    Lines = match.readlines()
    l_SRR={}
    matched_SRR = open('matching_SRR.txt', 'w')
    for line in Lines:
        matched_SRR.write(line.split(" ")[0] + "\n")
        line = line.strip()
        key = line.split(" ")[1] 
        l_SRR[key] = line.split(" ")[0]
    

    match.close()
    return l_SRR


def waitUntil(driver, by, xpath):
    timeout = 5
    try:
        element_present = EC.presence_of_element_located((by,xpath))
        return WebDriverWait(driver, timeout).until(element_present)
    except TimeoutError:
        print("PDT number not available")
    


def parse_cmdline_params(arg_list):
        """Parses commandline arguments.
        :param arg_list: Arguments to parse. Default is argv when called from the
        command-line.
        :type arg_list: list.
        """
        #Create instance of ArgumentParser
        options_parser = ArgumentParser(formatter_class=
                                        ArgumentDefaultsHelpFormatter)
        options_parser.add_argument('--SRR', dest='SRR', type=str, required=True,
                                    help="SRR number need to be searched")
        options_parser.add_argument('--folders', dest='folders', type=str, required=True,
                                    help="Destination folder to download")
        options_parser.add_argument('--file', dest='file', type=str, required=True,
                                    help="meta")
        opts = options_parser.parse_args(args=arg_list)
        return opts

    
if __name__ == "__main__":
    main(sys.argv)