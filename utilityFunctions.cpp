#include <sys/stat.h> 


std::string genRandomString(const int len) {
    char s[11];
    
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";

    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    s[len] = 0;
    std::string str(s);
    return( str );
}

bool FileExists(string strFilename) { 
  struct stat stFileInfo; 
  bool blnReturn; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if(intStat == 0) { 
    // We were able to get the file attributes 
    // so the file obviously exists. 
    cout << "Temp file already exists, trying again" << endl;
    blnReturn = true; 
  } else { 
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    blnReturn = false; 
  } 
  return(blnReturn); 
}

// std::string createTmpQualitiesFile (std::string temp, std::string inFile, int readLength, std::string qualfile) {
//     
//     bool done = false;
//     std::string tmpFile;
//     vector<string> quals;
//     static int fileNum = 1;
//     
//     if(qualfile == "") {
//         while(!done) {
//             tmpFile = temp + "/" + genRandomString(10) + ".tmp";
//             done = !FileExists(tmpFile);
//         }
//     }
//     else {
//         tmpFile = qualfile + boost::lexical_cast<std::string>(fileNum);
//     }
//         
//     std::ifstream f(inFile.c_str());      
//     if (!f.is_open())
//     {
//         std::cerr << "Unable to open qualities file: " << inFile << std::endl;
//         exit (EXIT_FAILURE); // error 
//     }
// 
//     std::ofstream outf(tmpFile.c_str());
//     if (!outf.is_open())
//     {
//         std::cerr << "Unable to create temporary file: " << tmpFile << std::endl;
//         exit (EXIT_FAILURE); // error 
//     }
//     
//     std::string buff;
//     std::cout << "Creating temporary qualities file: " << tmpFile << std::endl;
//     int count = 0;
// 
//     while (!f.eof()) {
//         count++;
//         getline(f,buff);
//         if(count == 16) {
//             outf << buff.substr((size_t) 0, (size_t) readLength) << std::endl;
//             count = 0;
//         }
//     }
//     f.close();
//     outf.close();
//     
//     return(tmpFile);
// }


std::string createTmpQualitiesFile2 (std::string temp, std::string inFile, int readLength, std::string qualfile) {
    
    bool done = false;
    std::string tmpFile;
    std::vector<std::string> quals;
    time_t start, end;
    static int fileNum = 1;
    
    if(qualfile == "") {
        while(!done) {
            tmpFile = temp + "/" + genRandomString(10) + ".tmp";
            done = !FileExists(tmpFile);
        }
    }
    else {
        tmpFile = temp + "/" + qualfile + boost::lexical_cast<std::string>(fileNum++) + ".tmp";
    }
        
    std::ifstream f(inFile.c_str());      
    if (!f.is_open())
    {
        std::cerr << "Unable to open qualities file: " << inFile << std::endl;
        exit (EXIT_FAILURE); // error 
    }

    std::string buff;
    std::cout << "Creating temporary qualities file: " << tmpFile << std::endl;
    int count = 0;

    time(&start);
    while (!f.eof()) {
        count++;
        getline(f,buff);
        if(count == 16) {
            quals.push_back( buff.substr((size_t) 0, (size_t) readLength) );
            count = 0;
        }
    }
    f.close();
    time(&end);
    printf("Reading time: %.0f seconds\n", difftime(end, start));
    
    
    
    std::ofstream outf(tmpFile.c_str());
    if (!outf.is_open())
    {
        std::cerr << "Unable to create temporary file: " << tmpFile << std::endl;
        exit (EXIT_FAILURE); // error 
    }
    
    
    time(&start);
    for(unsigned int i=0; i<quals.size(); i++) {
        outf << quals[i] << std::endl;
    }
    outf.close();
    time(&end);
    printf("Writing time: %.0f seconds\n", difftime(end, start));
    
    return(tmpFile);
}




std::string trim(std::string& s,const std::string& drop = " ") {
    std::string r=s.erase(s.find_last_not_of(drop)+1);
    return r.erase(0,r.find_first_not_of(drop));
}

int processLine(std::string& s) {
        if( s.compare(0,1,">") == 0 ) {
            return(1);
        }
        else {
            trim(s, "\n");
            return(0);               
        }
}


std::string reverseComplement(std::string read) {
    
    std::reverse(read.begin(), read.end());
    
    std::string::iterator It = read.begin();

    while ( It != read.end() )
    {
        switch (*It) {
            case 'A':
                *It = 'T';
                break;
            case 'C':
                *It = 'G';
                break;
            case 'G':
                *It = 'C';
                break;
            case 'T':
                *It = 'A';
                break;
        }
        *It++;
    }
    return(read);
}

std::string interpolateQualities(std::string quals, int desiredReadLength, int otherReadLength) {
    
    int i, n, idx;
    std::string newQuals(desiredReadLength, ' ');
    int readLength = desiredReadLength + otherReadLength;
    
    //cout << readLength << endl;
    
    std::vector<int> vals (readLength, -1);
        
    for(i = 0; i < (quals.length()); i++) {
        idx = (int) round( ((double) i / (double) (quals.length()-1) ) * (readLength - 1));
        vals[ idx ] = (int) quals[i];
    }
    
    for(i = 1; i < readLength - 1; i++) {
        if(vals[i] == -1) {
            n = 1;
            while(vals[i + n] == -1) n++;
            vals[i] = (int) round( (vals[i-1] + (vals[i+n] * n)) / (n+1) );
        }
    }    
    
    for(i = 0; i < desiredReadLength; i++) {
        newQuals[i] = (char) vals[i];
    }
    
    return newQuals;
}  
    

long wallTime(timeval start, timeval end) {
    
    long mtime, seconds, useconds;
 
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    
    mtime = (long int) ((seconds) * 1 + useconds/1000.0) + 0.5;
    
    return mtime;
}
