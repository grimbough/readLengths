#include "readLengths.h"

#include "classes.cpp"
#include "sampling.cpp"
#include "utilityFunctions.cpp"
//#include "classes_MP.cpp"

Genome readGenome(std::string inFile) {
    
    bool record = false;
    std::string buff, name, seq; 
       
    std::cout << "Reading genome file:" << inFile << std::endl;

    std::ifstream f(inFile.c_str());
    if (!f.is_open()) {
        cout << "Error reading genome file" << endl;
        exit(EXIT_FAILURE); // error 
    }

    Genome gen;

    while (!f.eof()) {
        getline(f,buff);
        if( processLine(buff) && record == true ) {

            gen.addChromosome(name, seq);
            name = "";
            seq = "";
            name += buff;
        }
        else {
            if(!record) {
                name += buff;
            }
            else {
                seq += buff;
            }
        }
        record = true;
    }
    f.close();
    gen.addChromosome(name, seq);
    
    return gen;    
}



int main(int argc, char *argv[]) {
    
    
    int readLength1, readLength2;
    std::string qualitiesFile1, qualitiesFile2, temp;
    std::string outFile, inFile;
    int numberOfReads, insertSizeMean, insertSizeSD;
    int qualityAdjust;
    unsigned seed;
    double snpProb;
    
    time_t start, end;
    
    /** BOOST Program Options **/
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("n", po::value<int>(), "number of read pairs to generate (default 50000)")
        ("mean", po::value<int>(), "mean insert size (default 200)")
        ("sd", po::value<int>(), "standard deviation of insert sizes (default 60)")
        ("l1", po::value<int>(), "read length of first read")
        ("l2", po::value<int>(), "read length of second read")
        ("q1", po::value<string>(), "qualities to be sampled for first read")
        ("q2", po::value<string>(), "qualities to be sampled for second read")
        ("genome", po::value<string>(), "path to reference genome reads will be sampled from")
        ("output", po::value<string>(), "")
        ("temp", po::value<string>(), "")
        ("seed", po::value<int>(), "")
        ("solexa", po::value<int>(), "")
        ("snp_prob", po::value<double>(), "")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return EXIT_FAILURE;
    }
    
    if (vm.count("n")) {
        numberOfReads = vm["n"].as<int>();
    } else {
        numberOfReads = 50000;
    }
    
    if (vm.count("mean")) {
        insertSizeMean = vm["mean"].as<int>();
    } else {
        insertSizeMean = 200;
    }
    
    if (vm.count("sd")) {
        insertSizeSD = vm["sd"].as<int>();
    } else {
        insertSizeSD = 60;
    }
   
    if (vm.count("l1")) {
        readLength1 = vm["l1"].as<int>();
    } else {
        readLength1 = 50;
    }
    
    if (vm.count("l2")) {
        readLength2 = vm["l2"].as<int>();
    } else {
        readLength2 = 50;
    }
    
    if (vm.count("q1")) {
        qualitiesFile1 = vm["q1"].as<string>();
    } else {
        cerr << "Qualities files for read one not provided" << endl;
        return EXIT_FAILURE;
    }
    
    if (vm.count("q2")) {
        qualitiesFile2 = vm["q2"].as<string>();
    } else {
        cerr << "Qualities files for read two not provided" << endl;
        return EXIT_FAILURE;
    }
    
    if (vm.count("genome")) {
        inFile = vm["genome"].as<string>();
    } else {
        cerr << "No genome file provided" << endl;
        return EXIT_FAILURE;
    }
    
    if (vm.count("output")) {
        outFile = vm["output"].as<string>();
    } else {
        cerr << "No output file specified" << endl;
        return EXIT_FAILURE;
    }
    
    if (vm.count("temp")) {
        temp = vm["temp"].as<string>();
    } else {
        temp = "/tmp";
    }
    
    if (vm.count("seed")) {
        seed = vm["seed"].as<int>();
    } else {
        seed = static_cast<unsigned> (std::time(0));
    }
    
    if (vm.count("solexa")) {
        qualityAdjust = 64;
    } else {
        qualityAdjust = 33;
    }
    
    if (vm.count("snp_prob")) {
        snpProb = vm["snp_prob"].as<double>();
    } else {
        snpProb = (double) 1 / 1000;
    }
    /** END BOOST Program Options **/
    
    // seed the random number generator
    // this is only used for the tmp file names, so don't use the seed
    srand( std::time(0) );
        
    time(&start);
    std::string tmpQuals1 = createTmpQualitiesFile(temp, qualitiesFile1, readLength1);
    std::string tmpQuals2 = createTmpQualitiesFile(temp, qualitiesFile2, readLength2);    
    time(&end);
    printf("Elapsed time: %.0f seconds\n", difftime(end, start));

    time(&start);
    Genome gen = readGenome(inFile);
    time(&end);
    printf("Elapsed time: %.0f seconds\n", difftime(end, start));

    time(&start);
    Reads reads = gen.generateReads(numberOfReads, readLength1, readLength2, insertSizeMean, insertSizeSD, seed);
    time(&end);
    printf("Elapsed time: %.0f seconds\n", difftime(end, start));
    
    time(&start);
    reads.generateQualities(tmpQuals1, tmpQuals2, seed);
    time(&end);
    printf("Elapsed time: %.0f seconds\n", difftime(end, start));
   // reads.generateQualities(tmpQuals2, 2);
    
    time(&start);
    #if defined (_OPENMP)
        reads.applyQualities_MP(seed, qualityAdjust, snpProb);
    #else
        reads.applyQualities(seed, qualityAdjust, snpProb);
    #endif
    time(&end);
    printf("Elapsed time: %.0f seconds\n", difftime(end, start));
     
    reads.writeFASTQ(outFile);

    remove(tmpQuals1.c_str());
    remove(tmpQuals2.c_str());

    return 0;

}

