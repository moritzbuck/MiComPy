import os
from tqdm import tqdm
import os.path
import subprocess
import shutil

class orthoMCL(object):
    orthoMCL_path = "/home/moritz/share/orthomclSoftware-v2.0.9/bin/"
    db_cfg = "/home/moritz/repos/Pyscratches/20150610_orthomcl_tools/orthmcl_tools/orthomcl.config.template"
    mysql_line = "mysql --defaults-file=~/sandboxes/msb_5_1_73/mysql.cnf -u orthomcl --protocol=TCP -P 5173"

    def __init__(self, out_dir, proteoms, name = "MCL"):
        self.out_dir = out_dir
        # list of prokka faa file
        self.proteoms = proteoms
        self.name = name
        self.out_mcl = out_dir+ "final_clusters.tsv"
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
   
    def full_pipe(self):
        self.empty_all_tables()
        print "Installing:"
  #      self.install_schemes()
        self.make_compliant_fastaa()
        self.make_good_fasta()
        self.do_the_blast()
        self.parse_blastout()
        self.load_into_db()
        self.make_pairs()
        self.dump_pairs()
        self.run_mcl()
        self.parse_mcl()

    def pre_blast(self):
        self.make_compliant_fastaa()
        self.make_good_fasta()


    def post_blast(self):
        self.start_server()
#        self.empty_all_tables()
        self.install_schemes()
        self.load_into_db()
        self.make_pairs()
        self.dump_pairs()
        self.stop_server()
        self.run_mcl()
        self.parse_mcl()
        self.stop_server()
        
    def start_server(self):
        print "Start MySQL server"
        start_cmd = " ".join([ '/home/moritz/sandboxes/msb_5_1_73/start', "--myisam_max_sort_file_size=6400000000", "--myisam_sort_buffer_size=7000000000","--read_buffer_size=250000000" ])
        os.system(start_cmd)

    def stop_server(self):
        print "Stopping MySQL server"
        os.system('/home/moritz/sandboxes/msb_5_1_73/stop')
        
        
    def install_schemes(self):
        print "installing orthoMCL DB-schemes"
        script = " ".join([self.orthoMCL_path + "orthomclInstallSchema", self.db_cfg])
        os.system(script)

    def make_compliant_fastaa(self):
        print "making compliant fastaas"
        if not os.path.exists(self.out_dir + "temp"):
            os.makedirs(self.out_dir + "temp")
        for f in tqdm(self.proteoms):
            script = " ".join([self.orthoMCL_path + "orthomclAdjustFasta", f.split("/")[-1].split(".")[0], f, str(1)])
            os.system(script)
            shutils.move(f.split("/")[-1].split(".")[0] + ".fasta",self.out_dir + "temp/" )

    def make_good_fasta(self, min_len = 10, max_stop_freq=20):
        print "making filtered fastaa"
        script = " ".join([self.orthoMCL_path + "orthomclFilterFasta", self.out_dir + "temp/", str(min_len), str(max_stop_freq) ])
        os.system(script)
        shutils.move("goodProteins.fasta", self.out_dir)

    def do_the_blast(self, processors = 16):
        print "running blast"
        os.system(" ".join([ "formatdb", "-i", self.out_dir + "goodProteins.fasta",  "-p", "T"]))
        blastall_cmd = ["blastall", "-p", "blastp", "-F", "'m S'", "-v", 100000, "-b", 100000, "-z", 0,  "-e", 1e-5, "-m", 8, "-a", processors, "-d", self.out_dir + "goodProteins.fasta", "-i", self.out_dir + "goodProteins.fasta", ">", self.out_dir + "blast_all_v_all.tsv"]
        blastall_cmd = [str(c) for c in blastall_cmd]
        os.system(" ".join(blastall_cmd))

    def parse_blastout(self):
        print "parsing the blastout"
        script = " ".join([self.orthoMCL_path + "orthomclBlastParser", self.out_dir + "blast_all_v_all.tsv", self.out_dir + "temp/",  " > ", self.out_dir + "blast_all_v_all.parsed.tsv"])
        os.system(script)

        
    def load_into_db(self):
        print "loading into the DB"
        script = " ".join([self.orthoMCL_path + "orthomclLoadBlast", self.db_cfg, self.out_dir + "blast_all_v_all.parsed.tsv"])
        os.system(script)

    def empty_all_tables(self):
        print "clean db"
        os.system('echo "TRUNCATE orthomcl.SimilarSequences;" | ' + self.mysql_line)
        os.system('echo "TRUNCATE orthomcl.CoOrtholog;" | ' + self.mysql_line)
        os.system('echo "TRUNCATE orthomcl.InParalog;" | ' + self.mysql_line)
        os.system('echo "TRUNCATE orthomcl.InParalog2Way;" | ' + self.mysql_line)
#        os.system('echo "TRUNCATE orthomcl.InterTaxonMatch;" | ' + self.mysql_line)
        os.system('echo "TRUNCATE orthomcl.Ortholog;" | ' + self.mysql_line)
        

    def make_pairs(self):
        print "making pairs with db-magik"
        script = " ".join([self.orthoMCL_path + "orthomclPairs", self.db_cfg, self.out_dir + "pairs.log", "cleanup=no"])
        os.system(script)

    def dump_pairs(self):
        if os.path.exists( "pairs/"):
            shutils.rmtree("pairs/")
        if os.path.exists(self.out_dir + "pairs/"):
            shutils.rmtree("-r", self.out_dir + "pairs/")
        print "dumping pairs"
        os.system(" ".join(self.orthoMCL_path + "orthomclDumpPairsFiles"))
        shutils.move("mclInput", self.out_dir)
        shutils.move( "pairs", self.out_dir)

    def run_mcl(self):
        print "Running MCL"
        os.system(" ".join([self.out_dir + "mclInput",  "--abc",  "-I", 1.5, "-o", self.out_dir + "mclOutput"]))

    def parse_mcl(self):
        print "parsing MCL"
        os.system(" ".join([self.orthoMCL_path + "orthomclMclToGroups", self.name + "_", "00000000", "<", self.out_dir + "mclOutput", ">", self.out_mcl]) )
