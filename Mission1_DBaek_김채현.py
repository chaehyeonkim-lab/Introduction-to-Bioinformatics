import time
import sys
import math

print(time.ctime())

# file processing
def file(file_name):
    
    handle = open(file_name, "r")
    handle.readline()
    fasta = handle.read().replace('\n','')
    handle.close()

    fasta = fasta.upper()
    
    return fasta
# end function

# class for sequencing counting
class seq_count:

    def __init__(self, Nmer):
        self.Nmer = int(Nmer)
        self.merDict = {}
    # count N-mers
    def counting(self, fasta):

        for i in range (len(fasta) - (self.Nmer-1)):
            nt = fasta[i:i + self.Nmer]
            
            if 'N' in nt:
                continue
            elif nt in self.merDict:
                self.merDict[nt] += 1
            else:
                self.merDict[nt] = 1
                
        return self.merDict
    # end function
    
    # calculate relative frequencies
    def ratio(self):
        self.total_count = 0
        self.ratio_sum = 0
        p_NT = {}

        for key in self.merDict:
            self.total_count += self.merDict[key]

        for key in self.merDict:
            freq = self.merDict[key] / self.total_count

            if freq < 0 or freq > 1:
                sys.exit("Probability is out of range")
            p_NT[key] = freq
            self.ratio_sum += freq

        return p_NT
    # end function
# end class

# calculate log2 fold change
def FC_Cal(twomer_p_NT, onemer_p_NT):
    cal_result = {}
    sum_expected = 0

    for key in twomer_p_NT:
        p_1stNT = key[0]
        p_2ndNT = key[1]
        
        expected_freq = onemer_p_NT[p_1stNT] * onemer_p_NT[p_2ndNT]
        if expected_freq < 0 or expected_freq > 1:
            sys.exit("Probability is out of range")

        sum_expected += expected_freq

        log2FC = math.log2(twomer_p_NT[key]/expected_freq)

        cal_result[key] = {"expected frequency": expected_freq, "log2FC": log2FC}

    return cal_result, sum_expected
# end function

# main function
def main(species_name, species, folder):
    onemer = seq_count(1)
    twomer = seq_count(2)

    for ch in species:
        file_name = '/Users/kimchaehyeon/Documents/seoul national univ./univ 4-2/생물학자를 위한 전산학 개론/'+folder+'/chr'+ch+'.fa'
        
        fasta = file(file_name) #데이터 가공 완료

        onemer_Dict = onemer.counting(fasta)
        twomer_Dict = twomer.counting(fasta)
    
    pVal_onemer = onemer.ratio() #1mer 비율 계산
    pVal_twomer = twomer.ratio() #2mer 비율 계산
        
    freq_cal, sum_expected = FC_Cal(pVal_twomer, pVal_onemer)

    print (species_name)
    for key in sorted(onemer_Dict):
        print(key, onemer_Dict[key], pVal_onemer[key], sep='\t')
    print ("Sum", onemer.total_count, onemer.ratio_sum, sep = '\t')
    for key in sorted(twomer_Dict):
        print (key, twomer_Dict[key], pVal_twomer[key], freq_cal[key]["expected frequency"], freq_cal[key]["log2FC"], sep = '\t')            
    print("Sum", twomer.total_count, twomer.ratio_sum, sum_expected, sep = '\t')
#end main function

C_elegans = ['I','II','III','IV','V','X']
D_melanogaster = ['2L','2R','3L','3R','X','Y']    
Chicken = [str(i) for i in range (1, 29)] + ['W','Z']
Human = [str(i) for i in range (1, 23)] + ['M','X','Y']

main("C_elegans", C_elegans, 'ce11ChrFiles')
main("D_melanogaster", D_melanogaster, 'dm6ChrFiles')
main("Chicken", Chicken, 'galGal6ChrFiles')
main("Human", Human, 'hg38ChrFilesIncludingChrM')

print(time.ctime())
