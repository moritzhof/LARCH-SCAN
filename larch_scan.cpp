/*#######################################################################################################

   			  GlobDiversity -- LARCH-SCAN Fragmentation Algorithm
 
    				Created: Feb 4, 2020 by:
				DLR: Moritz Travis Hof, Luca Spataro
				WUR: Rene Jochem 

#########################################################################################################*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

//##################################################################################################
//################################# Matrix namespace: float ########################################
//##################################################################################################

namespace matrix {
   typedef std::vector<float> Vector;
   typedef std::vector<Vector> Matrix;
}

//############################### Important Constants ##############################
//--------------- TO DO: read in these arguments are execution ---------------------   // Completed 
//---------------   alpha: 10km (0.299573) and 5km (0.599146)  ---------------------
//##################################################################################

//float percentage = 95.0;
//float alpha = 0.299573;
//float maxDistance = ((log(100.0) - log(100.0 - percentage)) / alpha) * 1000.0;

//##################################################################################################
//################### I/O functions to for .flt and .hdr file format ###############################
//##################################################################################################

matrix::Matrix readFLTFile(std::string& filename, int& nrow, int& ncol,
                           float& maxDistance, float& cellsize, int& shift) {
  matrix::Matrix data(nrow + 2 * shift,
                      std::vector<float>(ncol + 2 * shift, 0));
  std::ifstream infile(filename, std::ios::in | std::ifstream::binary);

  std::set<float> values;
  int i, row = 0, col = 0;
  float n;
  if (infile.is_open()) {
    i = 0;
    while (infile.read((char*)&n, sizeof(float))) {
      row = i / ncol;
      col = i % ncol;

      values.insert(n);

      data[row + shift][col + shift] = n;

      i = i + 1;
    }

    infile.close();
    std::cout
        << "Reading in values from .flt file and creating matrix completed"
        << std::endl;
  }
  return data;
}

//##################################################################################################
void writeFLTfile(matrix::Matrix& data, std::string fname) {
  std::string filename = fname+".flt";
  std::ofstream out_vec(filename, std::ios::out | std::ofstream::binary);
  float n;
  if (out_vec.is_open()) {
    for (std::vector<std::vector<float>>::iterator it = data.begin();
         it != data.end(); ++it) {
      for (std::vector<float>::iterator it2 = (*it).begin(); it2 != (*it).end();
           ++it2) {
        n = *it2;
        out_vec.write((char*)&n, sizeof(float));
      }
    }
  }
  out_vec.close();
}

//#################################################################################################
std::tuple<int, int, float, float, float> readHDRfile(std::string filename) {
  int col, row;
  float cornerX, cornerY, cellsize;
  std::fstream infile(filename);
  std::string b;
  if (!infile) {
    std::cerr << " Could not open file:" << filename << std::endl;
  }
  char buffer[90];

  infile >> buffer;
  infile >> buffer;
  b = buffer;
  col = std::stol(b);

  infile >> buffer;
  infile >> buffer;
  b = buffer;
  row = std::stol(b);

  infile >> buffer;
  infile >> buffer;
  b = buffer;
  cornerX = std::stod(b);

  infile >> buffer;
  infile >> buffer;
  b = buffer;
  cornerY = std::stod(b);

  infile >> buffer;
  infile >> buffer;
  b = buffer;
  cellsize = std::stod(b);

  //  infile >> buffer;
  infile.close();
  std::cout << "Reading in values from .hdr file and creating tuple completed"
            << std::endl;
  return std::make_tuple(col, row, cornerX, cornerY, cellsize);
}

//###############################################################################################
int writeHDRFile(std::string filenameinput, int rows, int cols, float& rowcorner, float& colcorner, float& cellsize){

  std::string filename = filenameinput+".hdr";
  std::ofstream infl;
  infl.open(filename.c_str());
  if (!infl) {
    throw("can not write file: " + filename);
  }
  infl << "ncols ";
  infl << cols;
  infl << std::endl;
  infl << "nrows ";
  infl << rows;
  infl << std::endl;
  infl << "xllcorner ";
  infl << rowcorner;
  infl << std::endl;
  infl << "yllcorner ";
  infl << colcorner;
  infl << std::endl;
  infl << "cellsize ";
  infl << cellsize;
  infl << std::endl;
  infl << "NODATA_value ";
  infl << -9999;
  infl << std::endl;
  infl << "byteorder     LSBFIRST";
  infl << std::endl;

  infl.close();
  return 1;
}

//#####################################################################################################
// If necessary to print and look at the information in the flt file read.
void printData(matrix::Matrix& data) {
  for (const auto& info : data) {
    for (const auto& value : info) {
      std::cout << value << ' ';
    }
    std::cout << '\n';
  }
}

//######################################################################################################
void printHelp(char* argv[]){

        std::cout<<"Programm: "<<argv[0]<<std::endl;

	std::cout<<"Parameters: "<<std::endl;
	std::cout<<"-population_map: Population Map Rasterfile"<<std::endl;
	std::cout<<"-scan_map: Scan Map Rasterfile"<<std::endl;


	std::cout<<"-alpha: Alpha (float), default: 0.299573 "<<std::endl;
	std::cout<<"-percentage: percentage [0-100], default 95"<<std::endl;


}

//#######################################################################################################
//################### Main Function for the LARCH SCAN algorithm. #######################################
//#######################################################################################################
 
matrix::Matrix doCell(matrix::Matrix& movingWindow, matrix::Matrix& flt,
                      int& orows, int& ocols, float maxDistance,
                      float ncellsize, int shift) {
  int n = flt.size();
  int m = flt[0].size();
  int p = movingWindow[0].size();
  matrix::Matrix floatGrid(orows + 2 * shift, std::vector<float>(ocols + 2 * shift, 0));


  
  std::cout << "flt size " << flt.size() << "   " << flt[0].size() << std::endl;
  std::cout << "moving window size  " << movingWindow.size() << "   "
            << movingWindow[0].size() << std::endl;

  for (int i = 0; i < orows; ++i) {
    for (int j = 0; j < ocols; ++j) {
      auto value = flt[i + shift][j + shift];
      if (value > 0) {
        floatGrid[i + shift][j + shift] += movingWindow[0][0];

        // go over middle row and column
        for (int k = 1; k < movingWindow[0].size(); ++k) {
          auto v = value * movingWindow[0][k];
          floatGrid[i + shift][j + shift + k] += v;
          floatGrid[i + shift][j + shift - k] += v;
          floatGrid[i + shift + k][j + shift] += v;
          floatGrid[i + shift - k][j + shift] += v;
        }

        for (int l = 1; l < movingWindow.size(); ++l) {
          for (int k = 1; k < movingWindow[l].size(); ++k) {
            // if(l = 0 && k == 0) continue;
            auto v = value * movingWindow[l][k];
            floatGrid[i + shift + l][j + shift + k] += v;
            floatGrid[i + shift - l][j + shift + k] += v;
            floatGrid[i + shift + l][j + shift - k] += v;
            floatGrid[i + shift - l][j + shift - k] += v;
          }
        }
      }
    }
  }
  return floatGrid;  // result matrix;
}


double ScaleMovWindow(matrix::Matrix& FMovingWindow) {
  double FSumWindow = 0.0;
  //! calulation three quadrants excluding 0
  for (int r = 1; r < FMovingWindow.size(); r++) {
    for (int c = 1; c < FMovingWindow[r].size(); c++) {
      FSumWindow += FMovingWindow[r][c];
    }
  }
  FSumWindow *= 4.0;
  //! calulation fourth quadrant including 0
  double sum = 0.0;
  for (int c = 1; c < FMovingWindow[0].size(); c++) {
    sum += FMovingWindow[0][c];
  }
  FSumWindow += sum * 4;
  FSumWindow += FMovingWindow[0][0];

  //! normalise moving window
  for (int r = 0; r < FMovingWindow.size(); r++) {
    for (int c = 0; c < FMovingWindow[r].size(); c++) {
      FMovingWindow[r][c] = FMovingWindow[r][c] / FSumWindow;
    }
  }

  return FSumWindow;
}

//#################################################################################################
// 		The Moving Window creates a spares matrics with a maximum radius
//			 (depending on alpha) around a non-zero value.
//#################################################################################################
matrix::Matrix MakeMovingWindow(float maxDistance, float cellsize, float alpha) {
  //! Get stamp matrix size
  float stampSize = 1 + maxDistance / cellsize;
  std::cout << "moving window stamp size: " << stampSize << std::endl;

  matrix::Matrix FMovingWindow;

  //! Create alpha distance matrix
  matrix::Vector p;
  int r = 0;
  double rd = double(r) * cellsize;
  while (rd <= maxDistance) {
    //! Create the row and calulate the square row distance
    double d = rd;
    double con;
    matrix::Vector dv;
    con = exp(-alpha * (d / 1000.0));
    dv.push_back(con);
    FMovingWindow.push_back(dv);
    double rd2 = rd * rd;

    //! Fill the row and calulate the square column distance.
    //! Add cells as long as the cell distance is
    //! within or equal to the max. distance.
    int c = 1;
    double cd = double(c) * cellsize;
    d = sqrt(rd2 + (cd * cd));
    con = exp(-alpha * (d / 1000.0));
    while (d <= maxDistance) {
      FMovingWindow[r].push_back(con);
      c++;
      cd = double(c) * cellsize;
      d = sqrt(rd2 + (cd * cd));
      con = exp(-alpha * (d / 1000.0));
    }
    r++;
    rd = double(r) * cellsize;
  }

  ScaleMovWindow(FMovingWindow);

  //printData(FMovingWindow);

  return FMovingWindow;
}

//#################################################################################################
// Main exeuction to lauch algorithm in the intended sequence
matrix::Matrix execute(matrix::Matrix& flt, int& ncol, int& nrow,
                       float& cellsize, float& alpha, float maxDistance,
                       int shift) {
  matrix::Matrix imovingWindow, result;
  std::cout << "computation commensing " << std::endl;

  imovingWindow = MakeMovingWindow(maxDistance, cellsize, alpha);
 
  std::cout << "moving window complete" << std::endl;

  result = doCell(imovingWindow, flt, nrow, ncol, maxDistance, cellsize, shift);

  std::cout << "computation completed" << std::endl;
  return result;
}

  
//###############################################################################
//             			  MAIN 
//###############################################################################

auto main(int argc, char* argv[]) -> int {
//int main( int argc, char* argv[] ){
   std::string population_map;
   std::string scan_map;
	
   float percentage = 95.0;
   float alpha = 1.0;
	
 //###############################################################################
 //################## Read in execution arguments ################################ 
	
  std::vector<std::string> arguments;

  if (argc <= 4) {
	std::cerr<< "No Parameters given. Exiting"<<std::endl;
	printHelp(argv);
	exit(1);
			 }
  for(int i = 1; i<argc; ++i){
	arguments.push_back(argv[i]);
  }
	
  population_map = std::move(arguments[0]);
  scan_map       = std::move(arguments[1]);
  alpha          = std::move(std::stod(arguments[2]));
  percentage     = std::move(std::stod(arguments[3]));

  //#############################################################################
  std::cout << population_map << "\t" << scan_map << std::endl; 
  // NOTE: smart to have the hdr file and flt file with the same name
  // read in .dhr file
  matrix::Matrix result;
  auto hdrFile = readHDRfile(population_map+".hdr");
 
  // assigned values hdr file.
  auto ncol 	 = std::move(std::get<0>(hdrFile));
  auto nrow 	 = std::move(std::get<1>(hdrFile));
  auto cornerX   = std::move(std::get<2>(hdrFile));
  auto cornerY   = std::move(std::get<3>(hdrFile));
  auto ncellsize = std::move(std::get<4>(hdrFile));
 
 	
  float maxDistance = ((log(100.0) - log(100.0 - percentage)) / alpha) * 1000.0;

  int shift = 1 + maxDistance / ncellsize;
  cornerX = cornerX - shift*ncellsize;
  cornerY = cornerY - shift*ncellsize;  
   
 
  // read in .flt file: dependent on .hdr file. 
  std::string filename = population_map +".flt";
 
  auto flt = readFLTFile(filename, nrow, ncol, maxDistance, ncellsize, shift);


  // Execution of primary algorithm: calls the Larch_Scan method
  result = execute(flt, ncol, nrow, ncellsize, alpha, maxDistance, shift);

  std::cout << result.size() << '\t' << result[0].size() << std::endl;
	
  std::cout << alpha << std::endl;
  writeFLTfile(result, scan_map);
  writeHDRFile(scan_map, result.size(), result[0].size(), cornerX, cornerY, ncellsize);
  std::cout << "Simulation fully completed" << std::endl;

  return 0;
}
