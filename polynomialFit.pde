int size = 5; //matrix size (n by n)
int degree = 3; //degree of polynomial regression

polyFit myFit = new polyFit();

void setup() {
  size(1000, 360);
  smooth();
  noStroke();
}

void draw() {
  int[] arrayX = new int[]{1,2,3,4,5};
  int[] arrayY = new int[]{100,200,500,1200,3500};
  int size = arrayX.length;
  int degree = 3;
  
  myFit.run(arrayX, arrayY, size, degree);
  double[] coeff = myFit.getCoeff();
  
  for(int i = 0; i <= degree; ++i){
    text((float)coeff[i],10,(i+1)*30);
  }
}

class polyFit{
  
  private int _size;
  private int _degree = degree;
  private int _numEq = _degree+1;
  private int[] _xVals = new int[_size]; 
  private int[] _yVals = new int[_size];
  private double[] xSigmaVals = new double[2*_degree+1];
  private double[] ySigmaVals = new double[_degree+1];
  private double[][] _augmentedMatrix = new double[_degree+1][_degree+2];
  private double[] _finalCoeff = new double[_degree+1];
  
  //stores the values of sigma(xi), sigma(xi^2),sigma(xi^3)...sigma(xi^2n)
  private void buildXSigmaVals(){
    for(int i=0; i < 2*_degree+1; ++i){
      xSigmaVals[i] = 0;
      for(int j = 0; j<_size; ++j){
        xSigmaVals[i] = xSigmaVals[i] + pow(_xVals[j],i);
      }
    }
  }
  
  void setX(){
    for(int i = 0; i < _size; ++i){
      _xVals[i] = i+1;
    }
    _yVals[0] = 100;
    _yVals[1] = 200;
    _yVals[2] = 500;
    _yVals[3] = 1200;
    _yVals[4] = 3500;
    
  }
  
  //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
  private void buildYSigmaVals(){
    for(int i = 0; i < _degree+1; ++i){
      ySigmaVals[i] = 0;
      for(int j = 0; j < _size; ++j){
        ySigmaVals[i] = ySigmaVals[i] + pow(_xVals[j],i)*_yVals[j];
      }
    }
  }
  
  //bulds the augmented matrix, stores corresponding at the right, except for last column
  //loads values of y as the last column of augmented matrix
  private void buildAugmentedMatrix(){
    for(int i = 0; i <= _degree; ++i){
      for(int j = 0; j <= _degree; ++j){
        _augmentedMatrix[i][j] = xSigmaVals[i+j];
      }
    }
    for(int i = 0; i <=_degree; ++i){
      _augmentedMatrix[i][_degree+1] = ySigmaVals[i];
    }
  }
  
  //gaussian elimination starts with pivotization procedure
  //puts each row in order
  private void pivotization(){
    for(int i = 0; i < _numEq; ++i){
      for(int k = i + 1; k < _numEq; ++k){
        if(_augmentedMatrix[i][i] < _augmentedMatrix[k][i]){
          for(int j = 0; j <= _numEq; ++j){
            swap(i,j,k,j);
          }
        }
      }
    }
  }
  
  //gaussian elimination
  private void gaussElim(){
    for(int i = 0; i < _numEq-1; ++i){
      for(int k = i + 1; k < _numEq; ++k){
        double coeff = _augmentedMatrix[k][i]/_augmentedMatrix[i][i];
        for(int j = 0; j <= _numEq; ++j){
          _augmentedMatrix[k][j] = _augmentedMatrix[k][j] - coeff*_augmentedMatrix[i][j];
        }
      }
    }
  }
  
  //back substritions
  private void backSub(){
    for(int i = _numEq - 1; i >= 0; --i){
      _finalCoeff[i] = _augmentedMatrix[i][_numEq];
      for(int j = 0; j < _numEq; ++j){
        if(j != i){
          _finalCoeff[i] = _finalCoeff[i] - _augmentedMatrix[i][j]*_finalCoeff[j];
        }
      }
      _finalCoeff[i] = _finalCoeff[i]/_augmentedMatrix[i][i];
    }
  }

  
  //helper function to swap values in matrices
  private void swap(int lhsX, int lhsY, int rhsX, int rhsY){
    double temp = _augmentedMatrix[lhsX][lhsY];
    _augmentedMatrix[lhsX][lhsY] = _augmentedMatrix[rhsX][rhsY];
    _augmentedMatrix[rhsX][rhsY] = temp;
  }
  
  double[] getCoeff(){
    return _finalCoeff;
  }
  
  void run(int[] x, int[] y, int size, int degree){
     _size = size;
     _degree = degree;
     _xVals = x;
     _yVals = y;
     buildXSigmaVals();
     buildYSigmaVals();
     buildAugmentedMatrix();
     pivotization();
     gaussElim();
     backSub();
  }
  
  
}
