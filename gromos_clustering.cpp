#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <math.h>
#include <ctime>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <iostream>

int main(int argc,char*argv[]){
{
 // read input parameters from readline
 // matrix file prefix
 std::string MATRIX_ = argv[1];
 // number of matrix files
 unsigned nfiles = atoi(argv[2]);
 // distance cutoff
 double CUTOFF_ = atof(argv[3]);
 // number of frames
 unsigned nitems = atoi(argv[4]);
 // percentage of traj
 double MAXC_ = 0.6;

 // 0) preliminary stuff
 // create map neighbors and weights
 std::map< unsigned, std::vector< unsigned > > neighbors;
 std::map< unsigned, double > weights;

 // add to neighbor list the frame itself
 for (unsigned i = 0; i < nitems; ++i) {
    neighbors[i].push_back(i);
    weights[i] = 1.0;
 }

 // stuff needed to read files
 std::ifstream rmsdfile;
 // line
 std::string line;
 // buffer string
 std::string buf;     
 // temporary vector of double 
 std::vector<double> line_split(3, 0.0);
 // atom indexes
 unsigned i0, i1;
 // Root Mean Square Deviation
 double rmsd;
 // counter
 unsigned i;
 
 // 1) cycle on the number of matrix files
 for(unsigned ifile=0; ifile<nfiles; ++ifile){
  // prepare full name
  std::string fullname = MATRIX_;
  // if more than one file name, add suffix
  if(nfiles > 1){
   // convert ifile to string
   std::stringstream ss;
   ss << ifile;
   fullname = fullname + "." + ss.str();
  }
  // Read RMSD file
  rmsdfile.open(fullname.c_str());
  if (rmsdfile.is_open()) {
   // read line by line
   while ( getline (rmsdfile,line) )
   { 
      // split line into string separated by a space
      std::stringstream ss(line);
      // reset counter
      i=0;
      while (ss >> buf) {line_split[i]=atof(buf.c_str()); ++i;}
      // convert
      i0 = static_cast<unsigned>(line_split[0]);
      i1 = static_cast<unsigned>(line_split[1]);
      rmsd = line_split[2];
      // add frames within cutoff
      if (rmsd < CUTOFF_) {
       neighbors[i0].push_back(i1);
       neighbors[i1].push_back(i0);
       weights[i0] += 1.0;
       weights[i1] += 1.0;
     }  
   }
   rmsdfile.close();
  }
  else std::cout << "Unable to open file" << std::endl;
  // end of files reading
 }
 
 // prepare list of clusters
 std::vector< std::vector<unsigned> > clusters;

 // start iterative procedure
 double maxweight = 1.0;
 while (maxweight > 0.0) {

    // find frame with maximum number of neighbors (weight) 
    maxweight = -1.0;
    int icenter = -1;
    // iterate on map
    for (std::map< unsigned, double >::iterator it = weights.begin(); it != weights.end(); ++it){
      if (it->second > maxweight) {
        maxweight = it->second;
        icenter = it->first;
      }
    }
    // no more clusters to find
    if (maxweight < 0.) {
      break;
    }

    // create the new cluster
    std::vector<unsigned> newcluster = neighbors[icenter];
    clusters.push_back(newcluster);

    // now remove from pool
    // Two different methods: the efficiency depends on cluster size
    if(newcluster.size() > MAXC_ * nitems){
     // this is more efficient with big clusters
     for (unsigned i = 0; i < newcluster.size(); ++i) {
       // remove entry from neighbors and weights maps
       neighbors.erase(newcluster[i]);
       weights.erase(newcluster[i]);
     }
     // cycle on cluster members
     for (unsigned i = 0; i < newcluster.size(); ++i) {
       for (std::map< unsigned, std::vector< unsigned > >::iterator it = neighbors.begin(); it!=neighbors.end(); ++it){
          // find element
          std::vector<unsigned>::iterator iit =
               find((it->second).begin(), (it->second).end(), newcluster[i]);
          if(iit != (it->second).end()){
           (it->second).erase(iit);
           weights[it->first] -= 1.0;
          }
       }
     }
    } else {
     // this is more efficient with small clusters
     for (unsigned i = 0; i < newcluster.size(); ++i) {
      // cycle on neighbors of i-th cluster element, excluded itself (j=0)
      for (unsigned j=1; j < neighbors[newcluster[i]].size(); ++j){
         // remove newcluster[i] from neighbor list
         unsigned index = neighbors[newcluster[i]][j];
         std::vector<unsigned>::iterator it =
              find(neighbors[index].begin(), neighbors[index].end(), newcluster[i]);
         neighbors[index].erase(it);
         weights[index] -= 1.0;
      }
      // remove entry from neighbors and weights maps
      neighbors.erase(newcluster[i]);
      weights.erase(newcluster[i]);
     }
    }
 // end of iterations
 }

 // print out stuff
 std::cout << "NUMBER OF CLUSTERS " << clusters.size() << std::endl;
 // print cluster statistics
 // open final log file
 FILE * log_final;
 log_final = fopen ("log.dat","w");
 for (unsigned i=0; i<clusters.size(); ++i){
  fprintf (log_final, "ID %10u  POPULATION %10lu  CENTER %10u\n",i, clusters[i].size(), clusters[i][0]);
 }
 fclose (log_final);

 // trajectory file 
 FILE * log_traj;
 log_traj = fopen ("trajectory.dat","w");
 // prepare list of assignments
 std::vector<unsigned> assign;
 for(unsigned i=0; i<nitems; ++i) assign.push_back(0);
 // cycle on clusters
 for(unsigned i=0; i<clusters.size(); ++i){
  for(unsigned j=0; j<clusters[i].size(); ++j){
    assign[clusters[i][j]]=i;
  }
 }
 // print out
 for(unsigned i=0; i<assign.size(); ++i) fprintf (log_traj, "%10u %10u\n",i,assign[i]);
 
 fclose (log_traj);

}
}
