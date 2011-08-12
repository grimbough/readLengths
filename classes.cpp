#include "readLengths.h"

class Chromosome {
        std::string name;
	std::string sequence;
public:
	Chromosome(string n, string s);
	void printName()  { cout << name << "\n"; };
	void printSeq() { cout << sequence << "\n"; };
	int length() { return(sequence.length()); };
	string returnSequence() { return(sequence); };
	string returnName() { return(name); };
};

class Reads {
	int readLength1;
        int readLength2;
	int insertSizeMean;
	int insertSizeSD;
	vector<string> end1;
	vector<string> end2;
        vector<string> quals1;
        vector<string> quals2;
	vector<string> chrom1;
	vector<string> chrom2;
	vector<string> pos1;
	vector<string> pos2;
        vector<bool> strand1;
public:
	Reads(int l1, int l2, int insM, int insSD) { readLength1 = l1; readLength2 = l2; insertSizeMean = insM; insertSizeSD = insSD;};
	void addReads(string, string, string, string, string, string, int);
        void generateQualities(std::string qualsFile1, std::string qualsFile2, int seed);
        void generateQualities(std::string qualsFile1, std::string qualsFile2, int seed, int readLength1, int readLength2);
        void applyQualities(int seed, int qualityAdjust, double snpProb);
        void applyQualities_MP(int seed, int qualityAdjust, double snpProb);
	void writeFASTA();
	void writeFASTA(string fileName);
        void writeFASTQ(string fileName);
};


class Genome {
	int numChrom;
	vector<Chromosome> chromVec;
public:
	Genome();
	void addChromosome(string, string);
	void printChromosomes();
	Reads generateReads(int, int, int, int, int, int);
        double* chromosomeWeights();
};




/** ******************************************* **/

Chromosome::Chromosome(string n, string s) {
	name = n;
	sequence = s;
}


/** ******************************************* **/

void Reads::addReads(string read1, string read2, string chr1, string chr2, string p1, string p2, int firstReadStrand) {
        
	end1.push_back(read1);
	end2.push_back(read2);
	chrom1.push_back(chr1);
	chrom2.push_back(chr2);
	pos1.push_back(p1);
	pos2.push_back(p2);
        strand1.push_back((bool) firstReadStrand);
}

void Reads::generateQualities(std::string tmpQuals1, std::string tmpQuals2, int seed) {
    
    int nlines, pos;
    unsigned i;
    int qualsLength1, qualsLength2;
    std::string buff; 
    
    std::cout << "Generating Qualities" << std::endl;
    
    std::ifstream qualsf1(tmpQuals1.c_str());
    std::ifstream qualsf2(tmpQuals2.c_str());
    if (!qualsf1.is_open() | !qualsf2.is_open())
    {
        std::cerr << "error reading temporary qualities file" << std::endl;
        exit (EXIT_FAILURE); // error 
    }
    /* find the number of lines */
    ///first read one line from each file and find their lengths
    getline(qualsf1, buff);
    qualsLength1 = buff.length();
    getline(qualsf2, buff);
    qualsLength2 = buff.length();
    
    /// now go to the end of the file
    qualsf1.seekg(0, std::ios::end);
    qualsf2.seekg(0, std::ios::end);
        
    if ( (int) (qualsf1.tellg() / (qualsLength1 + 1)) != (int) (qualsf2.tellg() / (qualsLength2 + 1)) ) {
        cerr << "Number of lines in qualities files do not match" << endl;
        exit(EXIT_FAILURE);
    }
    
    nlines = (int) (qualsf2.tellg() / (qualsLength2 + 1));
    
    for(i = 0; i < end1.size(); i++) {
        int line = SampleUniform(0, nlines - 1, seed);
        pos = (qualsLength1 + 1) * line;
        
        qualsf1.seekg(pos, std::ios::beg);
        getline(qualsf1, buff);
        if(buff.length() < readLength1) {
            buff = interpolateQualities(buff, readLength1);
        }
        quals1.push_back(buff);
        
        pos = (qualsLength2 + 1) * line;
        qualsf2.seekg(pos, std::ios::beg);
        getline(qualsf2, buff);
        if(buff.length() < readLength2) {
            buff = interpolateQualities(buff, readLength2);
        }
        quals2.push_back(buff);
    }
    //cout << "Made it" << endl;
    qualsf1.close();
    qualsf2.close();
 
}


void Reads::applyQualities(int seed, int qualityAdjust, double snpProb) {
 
    std::string bases = "ACGT";
    unsigned i, j;
    int q, baseChoice;
    double p;
    
    cout << "Applying qualities" << endl;
    
    for(i = 0; i < end1.size(); i++) {
        for(j = 0; j < end1[i].length(); j++) {
            
            // get the quality value
            q = (int) quals1[i][j] - qualityAdjust;
            // convert into a probaility of a miscall
            p = pow(10, (double) -q / 10);
            if(p != 1)
                p += snpProb;
            
            if(p > 1) {
                cout << "End1: " << quals1[i] << endl;
                cout << i << "\t" << j << "\t" << p << "\t" << quals1[i][j] << endl;
                exit(EXIT_FAILURE);
            }
                
            if(SampleBernoulli(p, seed)) {
                baseChoice = SampleUniform(0, 3, seed);
                end1[i][j] = bases[baseChoice];
            }
        }
    }
             
    for(i = 0; i < end2.size(); i++) {
        for(j = 0; j < end2[i].length(); j++) {        
            // repeat for the second read
            q = (int) quals2[i][j] - qualityAdjust;
            p = pow(10, (double) -q / 10);
            if(p != 1)
                p += snpProb;
            
            if(p > 1) {
                 cout << "End: " << quals2[i] << endl;
                 cout << i << "\t" << j << "\t" << p << "\t" << quals2[i][j] << endl;
                 exit(EXIT_FAILURE);
             }

            if(SampleBernoulli(p, seed)) {
                baseChoice = SampleUniform(0, 3, seed);
                end2[i][j] = bases[baseChoice];
            }
            
        }
    }
    
}

void Reads::writeFASTA() {
	unsigned i;
	for(i = 0; i < end1.size(); i++) {
		cout << chrom1[i] << ":" << pos1[i] << "/1\n" << end1[i] << "\n";
	}
	
	for(i = 0; i < end1.size(); i++) {
		cout << chrom2[i] << ":" << pos1[i] << "/2\n" << end2[i] << "\n";
	}
}

void Reads::writeFASTA(string fileName) {
    
        unsigned i;
        ofstream file1;
        ofstream file2;
        
        file1.open("test1.fa");    
        for(i = 0; i < end1.size(); i++) {
                file1 << chrom1[i] << ":" << pos1[i] << "/1\n" << end1[i] << "\n";
        }
        file2.close();
        
        file2.open("test2.fa");
        for(i = 0; i < end1.size(); i++) {
                file2 << chrom2[i] << ":" << pos2[i] << "/2\n" << end2[i] << "\n";
        }
        file2.close();
}

void Reads::writeFASTQ(string fileName) {
    
        unsigned i;
        ofstream file1;
        ofstream file2;
        
        std::string name1 = fileName + "_1.fq";
        std::string name2 = fileName + "_2.fq";
        
        file1.open(name1.c_str());    
        for(i = 0; i < end1.size(); i++) {
                file1 << "@" << chrom1[i].substr(1) << ":" << pos1[i] << ":" << pos2[i] << ":";
                
                if(strand1[i])
                    file1 << "+";
                else
                    file1 << "-";
                
                file1 << "/1" << std::endl;
                file1 << end1[i] << std::endl;
                file1 << "+" << std::endl;
                file1 << quals1[i] << std::endl;
        }
        file2.close();
        
        file2.open(name2.c_str());
        for(i = 0; i < end1.size(); i++) {
                file2 << "@" << chrom2[i].substr(1) << ":" << pos1[i] << ":" << pos2[i] << ":";
                
                if(strand1[i])
                    file2 << "+";
                else
                    file2 << "-";
                
                file2 << "/2" << std::endl;
                file2 << end2[i] << std::endl;
                file2 << "+" << std::endl;
                file2 << quals2[i] << std::endl;
        }
        file2.close();
}

/** ******************************************* **/

Genome::Genome() {
	numChrom = 0;
}

void Genome::addChromosome (string name, string sequence) {
	numChrom++;
	chromVec.push_back(Chromosome(name, sequence));
}

void Genome::printChromosomes() {
	
	int i;
	for(i = 0; i < numChrom; i++) {
		chromVec[i].printName();
//chromVec[i].printSeq();
	}
}

double* Genome::chromosomeWeights() {
    
        int i;
        vector<double> chromLengths;                 
        double *weights = (double *) calloc(sizeof(double), numChrom); 

        /* get all the chromosome lengths */
        for(i = 0; i < numChrom; i++) {
            chromLengths.push_back( (double) chromVec[i].length());
        }
        
        double max = *max_element(chromLengths.begin(), chromLengths.end());

        for(i = 0; i < numChrom; i++) {
            chromLengths[i] /= max ;
        }
        
        double sum_of_elems = 0;
        //for(std::vector<double>::iterator j=chromLengths.begin(); j!=chromLengths.end(); ++j)
        //    sum_of_elems += *j;
        for(i = 0; i < numChrom; i++) 
            sum_of_elems += chromLengths[i];
        
        //cout << sum_of_elems << endl;
        
        for(i = 0; i < numChrom; i++) {
            weights[i] = chromLengths[i] / sum_of_elems;
            
        }
        
        return(weights);
}


Reads Genome::generateReads(int totalReads, int readLength1, int readLength2, int insertSizeMean, int insertSizeSD, int seed) {
	
	int i, chrom, pos;
	string end1, end2, start1, start2;
	size_t found1, found2;
	Reads reads = Reads(readLength1, readLength2, insertSizeMean, insertSizeSD);
        vector<int> chromLengths;
        
        std::cout << "Generating reads" << std::endl;
        
        double *weights = chromosomeWeights();      
        //double weights[1] = {1.0};
        
        /* get all the chromosome lengths */
        for(i = 0; i < numChrom; i++) {
            chromLengths.push_back(chromVec[i].length());
        }
        
		
        for(chrom = 0; chrom < numChrom; chrom++) {
                       
            int numReads = (int) floor( ( weights[chrom] * totalReads ) + 0.5);
            
            //cout << chrom << "\t" << numReads << endl;
        
            for(i = 0; i < numReads; i++) {
                    
                int insertSize = SampleNormal((double) insertSizeMean, (double) insertSizeSD, seed);
                    
                    /// choose whether the first read is from the +ve or -ve strand
                    int firstReadStrand = (int) SampleBernoulli(0.5, seed);
                    
                    if(firstReadStrand) { // if 1 first read is on positive strand
                        pos = SampleUniform(0, chromVec[chrom].length() - (readLength1 + readLength2 + insertSize + 1), seed);
                        end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                        end2 = chromVec[chrom].returnSequence().substr(pos + readLength1 + insertSize, readLength2);
                        end2 = reverseComplement(end2);
                    }
                    else { // if 0 first read -ve strand
                        pos = SampleUniform(readLength2 + insertSize, chromVec[chrom].length() - (readLength1 + 1), seed);
                        end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                        end1 = reverseComplement(end1);
                        end2 = chromVec[chrom].returnSequence().substr(pos - insertSize - readLength2, readLength2);
                    }
                    
                    found1 = end1.find("N");
                    found2 = end2.find("N");
                    
                    while(found1 != string::npos || found2 != string::npos) {
                            //pos = rand() % (chromVec[chrom].length()  - (readLength1 + readLength2 + insertSize));
                            //pos = SampleUniform(0,  chromVec[chrom].length() - (readLength1 + readLength2 + insertSize + 1), seed);
                            
                            //end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                            //end2 = chromVec[chrom].returnSequence().substr(pos + readLength1 + insertSize, readLength2);
                            
                        if(firstReadStrand) { // if 1 first read is on positive strand
                            pos = SampleUniform(0, chromVec[chrom].length() - (readLength1 + readLength2 + insertSize + 1), seed);
                            end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                            end2 = chromVec[chrom].returnSequence().substr(pos + readLength1 + insertSize, readLength2);
                            end2 = reverseComplement(end2);
                        }
                        else { // if 0 first read -ve strand
                            pos = SampleUniform(readLength2 + insertSize, chromVec[chrom].length() - (readLength1 + 1), seed);
                            end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                            end1 = reverseComplement(end1);
                            end2 = chromVec[chrom].returnSequence().substr(pos - insertSize - readLength2, readLength2);
                        }
                        
                        found1 = end1.find("N");
                        found2 = end2.find("N");
                    }
                    
                    start1 = boost::lexical_cast< string >( pos + 1 );
                    if(firstReadStrand)
                        start2 = boost::lexical_cast< string >( pos + readLength1 + insertSize + 1 );
                    else
                        start2 = boost::lexical_cast< string >( pos - readLength2 - insertSize + 1 );
                    
                    reads.addReads(end1, end2, chromVec[chrom].returnName(), chromVec[chrom].returnName(), start1, start2, firstReadStrand);
                    
            }
            
        }
        
        free(weights);

	return(reads);
}
