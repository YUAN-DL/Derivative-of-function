  #include<iostream>
  #include<complex>
  #include<fftw3.h>
  #include<cstdlib>
  #include<cmath>
  #include<fstream>
  #define PI acos(-1)
  #define N 2048//number of points
  using namespace std;
  struct point
  {
      double x;
      double y;
  };
  class Function
  {
  private:
      point *points_arr;
      int number_points;
      double L;
      double h;
  public:
      Function(point *arg_p,int len_p,double L,double h);
      double* first_Derivative_Finite_difference_method();
      double* second_Derivative_Finite_difference_method();
      double* first_Derivative_Fourier_spectrum_method();
      double* second_Derivative_Fourier_spectrum_method();
      ~Function();
  };
  Function::Function(point *arg_p,int len_p,double arg_L,double arg_h)
  {
      if(len_p>0)
      { 
        number_points=len_p;
        L=arg_L;
        h=arg_h;
        points_arr=new point[len_p];
        for (int i = 0; i < len_p; i++)
        {
            points_arr[i].x=arg_p[i].x;
            points_arr[i].y=arg_p[i].y;
        }
      }
  }
  double* Function::first_Derivative_Finite_difference_method()
  {
      double *first_Derivative=new double[number_points];
      first_Derivative[0]=(points_arr[1].y-points_arr[0].y)/h;
      for (int i =1;i<number_points-1;i++)
      {
          first_Derivative[i]=(points_arr[i+1].y-points_arr[i-1].y)/(2*h);
      }
      first_Derivative[number_points-1]=(points_arr[number_points-1].y-points_arr[number_points-2].y)/h;
      return  first_Derivative;
  }
  double* Function::second_Derivative_Finite_difference_method()
  {
      double *second_Derivative =new double[number_points];
      double x=h*h;
      second_Derivative[0]=0.;
      for (int i = 1; i < number_points-1; i++)
      {     
          second_Derivative[i]=(points_arr[i-1].y-2*points_arr[i].y+points_arr[i+1].y)/x;
      }
      return second_Derivative;
  }
  double* Function::first_Derivative_Fourier_spectrum_method()
  {
      double *in;
      fftw_complex *out;
      fftw_plan p;
      in=(double *)fftw_malloc(sizeof(double)*N);
      out=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
      for(int i=0; i<N; i++)
      {
        in[i]=points_arr[i].y;
      }
      p=fftw_plan_dft_r2c_1d(N,in,out,FFTW_ESTIMATE);
      fftw_execute(p);
      double *k =new double[N];
      k[0]=0.0;
      for (int i = 1; i < N/2; i++)
      {
        k[i]=k[i-1]+(2*PI/L);    
       // k[i+N/2]=k[i]-N/2*(2*PI/L); 
      }
   // k[N/2]=k[0]-(N/2)*(2*PI/L);
      for (int i =0; i < N/2; i++)
      {   
        double a=-k[i]*out[i][1];
        double b=k[i]*out[i][0];
        out[i][0]=a;
        out[i][1]=b;
      }
      double *ifft;
      ifft=(double *)fftw_malloc(sizeof(double)*N);
      p=fftw_plan_dft_c2r_1d(N,out,ifft,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
      fftw_free(in);
      fftw_free(out);
      return ifft;
  }
  double* Function::second_Derivative_Fourier_spectrum_method()
  {
      double *in;
      fftw_complex *out;
      fftw_plan p;
      in=(double *)fftw_malloc(sizeof(double)*N);
      out=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
      for(int i=0; i<N; i++)
      {
        in[i]=points_arr[i].y;
      }
      p=fftw_plan_dft_r2c_1d(N,in,out,FFTW_ESTIMATE);
      fftw_execute(p);
      double *k =new double[N];
      k[0]=0.0;
      for (int i = 1; i < N/2; i++)
      {
        k[i]=k[i-1]+(2*PI/L);    
      }
      for (int i =0; i < N/2; i++)
      {   
        double a=-k[i]*k[i]*out[i][0];
        double b=-k[i]*k[i]*out[i][1];
        out[i][0]=a;
        out[i][1]=b;
      }
      double *ifft;
      ifft=(double *)fftw_malloc(sizeof(double)*N);
      p=fftw_plan_dft_c2r_1d(N,out,ifft,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
      fftw_free(in);
      fftw_free(out);
      return ifft;
  }
  Function::~Function()
  {
      delete[] points_arr;
  }
  void test1()
  {
      cout<<"Test1 first_Derivative"<<endl;
      point *p1=new point[N];
      int len=N;
      double L=2*PI;
      double h=L/N;
      p1[0].x=0;
      p1[0].y=exp(sin(p1[0].x));
      for (int i = 1; i < N; i++)
      {
          p1[i].x=p1[i-1].x+h;
          p1[i].y=exp(sin(p1[i].x));
      }
      Function F1(p1,N,L,h);
      double *first_Derivative_Finite;
      first_Derivative_Finite=F1.first_Derivative_Finite_difference_method();
      double *Analytical_solution_first=new double[N];
      for (int i = 0; i < N; i++)
      {
          Analytical_solution_first[i]=cos(i*h)*exp(sin(i*h));
      }
      cout<<"Analytical solution - first_Derivative："<<endl;
      for (int i = 0; i < N; i++)
      {
        // cout<<Analytical_solution_first[i]<<" ";
      }
      cout<<endl;
      cout<<endl;
      cout<<endl;
      cout<<"Finite difference solution - first_Derivative："<<endl;
      ofstream OutFile;
      OutFile.open("Test1_1_Finite_difference.txt");
      for (int i = 0; i < N; i++)
      {     
          
         // cout<<first_Derivative_Finite[i]<<" ";
          OutFile<<p1[i].x<<" "<<first_Derivative_Finite[i]<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *spectra_Methods_1=F1.first_Derivative_Fourier_spectrum_method();
      cout<<"Fourier spectra solution - first_Derivative："<<endl;
      OutFile.open("Test1_1_Fourier_spectra.txt");
      for (int i = 0; i < N; i++)
      {
          //cout<<spectra_Methods_1[i]/N<<" ";
          OutFile<<p1[i].x<<" "<<spectra_Methods_1[i]/N<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<"Test1 second_Derivative"<<endl;
      
      double *Analytical_solution_second=new double[N];
      for (int i = 0; i < N; i++)
      {
          Analytical_solution_second[i]=-sin(i*h)*exp(sin(i*h))+cos(i*h)*cos(i*h)*exp(sin(i*h));
      }
      cout<<"Analytical solution - second_Derivative："<<endl;
      for (int i = 0; i < N; i++)
      {
        // cout<<Analytical_solution_second[i]<<" ";
      }
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *error_first_Derivative_Finite=new double[N];
      for (int i = 0; i < N; i++)
      {
          error_first_Derivative_Finite[i]=abs(first_Derivative_Finite[i]-Analytical_solution_first[i]);
      }
      double error_max=-1.0;
      for (int i = 0; i < N; i++)
      {
          if(error_first_Derivative_Finite[i]>error_max)
            error_max=error_first_Derivative_Finite[i];
      }
      cout<<"The max error of Finite difference solution - first_Derivative: "<<error_max<<endl;
      delete[] error_first_Derivative_Finite;
      double *error_FFT=new double[N];
      for (int i = 0; i < N; i++)
      {
          error_FFT[i]=abs(spectra_Methods_1[i]/N-Analytical_solution_first[i]);
      }
      error_max=-1.0;
      for (int i = 0; i < N; i++)
      {
          if(error_FFT[i]>error_max)
            error_max=error_FFT[i];
      }
      cout<<"The max error of Fourier spectra solution - first_Derivative: "<<error_max<<endl;
      delete[]error_FFT;
    
      double *second_Derivative_Finite=F1.second_Derivative_Finite_difference_method();
      cout<<"Finite difference solution - second_Derivative："<<endl;

      OutFile.open("Test1_2_Finite_difference.txt");
      for (int i = 0; i < N-1; i++)
      {
          //cout<<second_Derivative_Finite[i]<<" ";
          OutFile<<p1[i].x<<" "<<second_Derivative_Finite[i]<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *spectra_Methods_2=F1.second_Derivative_Fourier_spectrum_method();
      
      cout<<"Fourier spectra solution - second_Derivative："<<endl;

      OutFile.open("Test1_2_Fourier_spectra.txt");
      for (int i = 1; i < N-1; i++)
      {
         //cout<<spectra_Methods_2[i]/N<<" ";
         OutFile<<p1[i].x<<" "<<spectra_Methods_2[i]/N<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *error_second_Derivative_Finite=new double[N];
      for (int i = 1; i < N-1; i++)
      {
          error_second_Derivative_Finite[i]=abs(second_Derivative_Finite[i]-Analytical_solution_second[i]);
      }
      error_max=-1.0;
      for (int i = 0; i < N-1; i++)
      {
          if(error_second_Derivative_Finite[i]>error_max)
            error_max=error_second_Derivative_Finite[i];
      }
      cout<<"The max error of Finite difference solution - second_Derivative: "<<error_max<<endl;
      delete[] error_second_Derivative_Finite;
      double *error_FFT_2=new double[N];
      for (int i = 1; i < N-1; i++)
      {
          error_FFT_2[i]=abs(spectra_Methods_2[i]/N-Analytical_solution_second[i]);
      }
      error_max=-1.0;
      for (int i = 1; i < N-1; i++)
      {
          if(error_FFT_2[i]>error_max)
            error_max=error_FFT_2[i];
      }
      cout<<"The max error of Fourier spectra solution - second_Derivative: "<<error_max<<endl;
      delete[]error_FFT_2;


      fftw_free(spectra_Methods_1);    
      fftw_free(spectra_Methods_2);
      delete[] p1;
      delete[] first_Derivative_Finite;
      delete[] Analytical_solution_first;
      delete[] second_Derivative_Finite;
      delete[] Analytical_solution_second;
  }
  void test2()
  {
      cout<<"Test2 first_Derivative"<<endl;
      point *p1=new point[N];
      int len=N;
      double L=2*PI;
      double h=L/N;
      p1[0].x=0;
      p1[0].y=exp(cos(p1[0].x));
      for (int i = 1; i < N; i++)
      {
          p1[i].x=p1[i-1].x+h;
          p1[i].y=exp(cos(p1[i].x));
      }
      Function F1(p1,N,L,h);
      double *first_Derivative_Finite;
      first_Derivative_Finite=F1.first_Derivative_Finite_difference_method();
      double *Analytical_solution_first=new double[N];
      for (int i = 0; i < N; i++)
      {
          Analytical_solution_first[i]=-sin(i*h)*exp(cos(i*h));
      }
      cout<<"Analytical solution - first_Derivative："<<endl;
      for (int i = 0; i < N; i++)
      {
         //cout<<Analytical_solution_first[i]<<" ";
      }
      cout<<endl;
      cout<<endl;
      cout<<endl;
      cout<<"Finite difference solution - first_Derivative："<<endl;
      ofstream OutFile;
      OutFile.open("Test2_1_Finite_difference.txt");
      for (int i = 0; i < N; i++)
      {     
         // cout<<first_Derivative_Finite[i]<<" ";
          OutFile<<p1[i].x<<" "<<first_Derivative_Finite[i]<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *spectra_Methods_1=F1.first_Derivative_Fourier_spectrum_method();
      cout<<"Fourier spectra solution - first_Derivative："<<endl;
      OutFile.open("Test2_1_Fourier_spectra.txt");
      for (int i = 0; i < N; i++)
      {
         // cout<<spectra_Methods_1[i]/N<<" ";
          OutFile<<p1[i].x<<" "<<spectra_Methods_1[i]/N<<endl;
      }
	  OutFile.close();
      cout<<endl;
      double *error_first_Derivative_Finite=new double[N];
      for (int i = 0; i < N; i++)
      {
          error_first_Derivative_Finite[i]=abs(first_Derivative_Finite[i]-Analytical_solution_first[i]);
      }
      double error_max=-1.0;
      for (int i = 0; i < N; i++)
      {
          if(error_first_Derivative_Finite[i]>error_max)
            error_max=error_first_Derivative_Finite[i];
      }
      cout<<"The max error of Finite difference solution - first_Derivative: "<<error_max<<endl;
      delete[] error_first_Derivative_Finite;
      double *error_FFT=new double[N];
      for (int i = 0; i < N; i++)
      {
          error_FFT[i]=abs(spectra_Methods_1[i]/N-Analytical_solution_first[i]);
      }
      error_max=-1.0;
      for (int i = 0; i < N; i++)
      {
          if(error_FFT[i]>error_max)
            error_max=error_FFT[i];
      }
      cout<<"The max error of Fourier spectra solution - first_Derivative: "<<error_max<<endl;
      delete[] error_FFT;
      cout<<"Test2 second_Derivative"<<endl;
      
      double *Analytical_solution_second=new double[N];
      for (int i = 0; i < N; i++)
      {
          Analytical_solution_second[i]=-cos(i*h)*exp(cos(i*h))+sin(i*h)*sin(i*h)*exp(cos(i*h));
      }
      cout<<"Analytical solution - second_Derivative："<<endl;
      for (int i = 0; i < N; i++)
      {
        // cout<<Analytical_solution_second[i]<<" ";
      }
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *second_Derivative_Finite=F1.second_Derivative_Finite_difference_method();
      cout<<"Finite difference solution - second_Derivative："<<endl;
      OutFile.open("Test2_2_Finite_difference.txt");
      for (int i = 0; i < N-1; i++)
      {
         // cout<<second_Derivative_Finite[i]<<" ";
          OutFile<<p1[i].x<<" "<<second_Derivative_Finite[i]<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *spectra_Methods_2=F1.second_Derivative_Fourier_spectrum_method();

      cout<<"Fourier spectra solution - second_Derivative："<<endl;
      OutFile.open("Test2_2_Fourier_spectra.txt");
      for (int i = 0; i < N-1; i++)
      {
         //cout<<spectra_Methods_2[i]/N<<" ";
         OutFile<<p1[i].x<<" "<<spectra_Methods_2[i]/N<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *error_second_Derivative_Finite=new double[N];
      for (int i = 1; i < N-1; i++)
      {
          error_second_Derivative_Finite[i]=abs(second_Derivative_Finite[i]-Analytical_solution_second[i]);
      }
      error_max=-1.0;
      for (int i = 1; i < N-1; i++)
      {
          if(error_second_Derivative_Finite[i]>error_max)
            error_max=error_second_Derivative_Finite[i];
      }
      cout<<"The max error of Finite difference solution - second_Derivative: "<<error_max<<endl;
      delete[] error_second_Derivative_Finite;
      double *error_FFT_2=new double[N];
      for (int i = 1; i < N-1; i++)
      {
          error_FFT_2[i]=abs(spectra_Methods_2[i]/N-Analytical_solution_second[i]);
      }
      error_max=-1.0;
      for (int i = 1; i < N-1; i++)
      {
          if(error_FFT_2[i]>error_max)
          {
            error_max=error_FFT_2[i];
          }
      }
      cout<<"The max error of Fourier spectra solution - second_Derivative: "<<error_max<<endl;
      delete[]error_FFT_2;


      fftw_free(spectra_Methods_1);    
      fftw_free(spectra_Methods_2);    
      delete[] p1;
      delete[] first_Derivative_Finite;
      delete[] Analytical_solution_first;
      delete[] second_Derivative_Finite;
      delete[] Analytical_solution_second;
  }  
  void test3()
  {
      cout<<"Test3 first_Derivative "<<endl;
      point *p1=new point[N];
      int len=N;
      double L=2*PI;
      double h=L/N;
      p1[0].x=0;
      p1[0].y=5*(sin(p1[0].x)+PI);
      for (int i = 1; i < N; i++)
      {
          p1[i].x=p1[i-1].x+h;
          p1[i].y=5*(sin(p1[i].x)+PI);
      }
      Function F1(p1,N,L,h);
      double *first_Derivative_Finite;
      first_Derivative_Finite=F1.first_Derivative_Finite_difference_method();
      double *Analytical_solution_first=new double[N];
      for (int i = 0; i < N; i++)
      {
          Analytical_solution_first[i]=5*(cos(p1[i].x));
      }
      cout<<"Analytical solution -first_Derivative："<<endl;
      for (int i = 0; i < N; i++)
      {
        // cout<<Analytical_solution_first[i]<<" ";
      }
      cout<<endl;
      cout<<endl;
      cout<<endl;
      cout<<"Finite difference solution -first_Derivative："<<endl;
      ofstream OutFile;
      OutFile.open("Test3_1_Finite_difference.txt");
      for (int i = 0; i < N; i++)
      {     
         // cout<<first_Derivative_Finite[i]<<" ";
          OutFile<<p1[i].x<<" "<<first_Derivative_Finite[i]<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *spectra_Methods_1=F1.first_Derivative_Fourier_spectrum_method();
      cout<<"Fourier spectra solution -first_Derivative："<<endl;
      OutFile.open("Test3_1_Fourier_spectra.txt");
      for (int i = 0; i < N; i++)
      {
         // cout<<spectra_Methods_1[i]/N<<" ";
          OutFile<<p1[i].x<<" "<<spectra_Methods_1[i]/N<<endl;
      }
	  OutFile.close();
      cout<<endl;
      double *error_first_Derivative_Finite=new double[N];
      for (int i = 0; i < N; i++)
      {
          error_first_Derivative_Finite[i]=abs(first_Derivative_Finite[i]-Analytical_solution_first[i]);
      }
      double error_max=-1.0;
      for (int i = 0; i < N; i++)
      {
          if(error_first_Derivative_Finite[i]>error_max)
            error_max=error_first_Derivative_Finite[i];
      }
      cout<<"The max error of Finite difference solution - first_Derivative: "<<error_max<<endl;
      delete[] error_first_Derivative_Finite;
      double *error_FFT=new double[N];
      for (int i = 0; i < N; i++)
      {
          error_FFT[i]=abs(spectra_Methods_1[i]/N-Analytical_solution_first[i]);
      }
      error_max=-1.0;
      for (int i = 0; i < N; i++)
      {
          if(error_FFT[i]>error_max)
            error_max=error_FFT[i];
      }
      cout<<"The max error of Fourier spectra solution - first_Derivative: "<<error_max<<endl;
      delete[] error_FFT;

      cout<<"Test3 second_Derivative"<<endl;
      
      double *Analytical_solution_second=new double[N];
      for (int i = 0; i < N; i++)
      {
          Analytical_solution_second[i]=-5*sin(i*h);
      }
      cout<<"Analytical solution - second_Derivative："<<endl;
      for (int i = 0; i < N; i++)
      {
        // cout<<Analytical_solution_second[i]<<" ";
      }
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *second_Derivative_Finite=F1.second_Derivative_Finite_difference_method();
      cout<<"Finite difference solution - second_Derivative："<<endl;
      OutFile.open("Test3_2_Finite_difference.txt");
      for (int i = 0; i < N-1; i++)
      {
        //  cout<<second_Derivative_Finite[i]<<" ";
          OutFile<<p1[i].x<<" "<<second_Derivative_Finite[i]<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *spectra_Methods_2=F1.second_Derivative_Fourier_spectrum_method();

      cout<<"Fourier spectra solution - second_Derivative："<<endl;
      OutFile.open("Test3_2_Fourier_spectra.txt");
      for (int i = 0; i < N; i++)
      {
        // cout<<spectra_Methods_2[i]/N<<" ";
         OutFile<<p1[i].x<<" "<<spectra_Methods_2[i]/N<<endl;
      }
	  OutFile.close();
      cout<<endl;
      cout<<endl;
      cout<<endl;
      double *error_second_Derivative_Finite=new double[N];
      for (int i = 1; i < N; i++)
      {
          error_second_Derivative_Finite[i]=abs(second_Derivative_Finite[i]-Analytical_solution_second[i]);
      }
      error_max=-1.0;
      for (int i = 1; i < N-1; i++)
      {
          if(error_second_Derivative_Finite[i]>error_max)
            error_max=error_second_Derivative_Finite[i];
      }
      cout<<"The max error of Finite difference solution - second_Derivative: "<<error_max<<endl;
      delete[] error_second_Derivative_Finite;
      double *error_FFT_2=new double[N];
      for (int i = 1; i < N-1; i++)
      {
          error_FFT_2[i]=abs(spectra_Methods_2[i]/N-Analytical_solution_second[i]);
      }
      error_max=-1.0;
      for (int i = 1; i < N-1; i++)
      {
          if(error_FFT_2[i]>error_max)
            error_max=error_FFT_2[i];
      }
      cout<<"The max error of Fourier spectra solution - second_Derivative: "<<error_max<<endl;
      delete[]error_FFT_2;      
      fftw_free(spectra_Methods_1);    
      fftw_free(spectra_Methods_2);    
      delete[] p1;
      delete[] first_Derivative_Finite;
      delete[] Analytical_solution_first;
      delete[] second_Derivative_Finite;
      delete[] Analytical_solution_second;
  }  
  int main(int argc, char** argv)
  {
      test1();
      test2();
      test3();
      return 0;
  }