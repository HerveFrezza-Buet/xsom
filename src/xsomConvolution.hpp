#pragma once

#include <fftw3.h>

namespace xsom {
  
  namespace tab {
    namespace fft {
      class Convolution {
      protected:
	typedef struct Workspace
	{
	  double * in_src, *out_src, *in_kernel, *out_kernel, *src_fft;
	  unsigned int w_fftw, h_fftw;
	  double * dst_fft;
	  double * dst; // The array containing the result
	  fftw_plan p_forw_src;
	  fftw_plan p_forw_kernel;
	  fftw_plan p_back;

	} Workspace;


	unsigned int width,height;
	unsigned int size;
	std::vector<double> data;
      public:
	mutable Workspace ws;
      protected:
	unsigned int kernel_width,kernel_height; 
	double* kernel;
	std::array<unsigned int,7> FFTW_FACTORS;

	void fftw_compute_kernel() {
	  double * ptr, *ptr_end;
	  for(ptr = ws.in_kernel, ptr_end = ws.in_kernel + ws.h_fftw*ws.w_fftw ; ptr != ptr_end ; ++ptr)
	    *ptr = 0.0;

	  for(unsigned int i = 0 ; i < kernel_height ; ++i)
	    for(unsigned int j = 0 ; j < kernel_width ; ++j)
	      ws.in_kernel[(i%ws.h_fftw)*ws.w_fftw+(j%ws.w_fftw)] += kernel[i*kernel_width + j];

	  fftw_execute(ws.p_forw_kernel);
	}

	void fftw_circular_convolution() const {
	  double * ptr, *ptr_end, *ptr2, *ptr3;

	  // Reset the content of ws.in_src
	  for(ptr = ws.in_src, ptr_end = ws.in_src + ws.h_fftw*ws.w_fftw ; ptr != ptr_end ; ++ptr)
	    *ptr = 0.0;

	  // Then we build our periodic signals
	  for(unsigned int i = 0 ; i < height ; ++i)
	    for(unsigned int j = 0 ; j < width ; ++j)
	      ws.in_src[(i%ws.h_fftw)*ws.w_fftw+(j%ws.w_fftw)] += data[i*width + j];

	  // And we compute their packed FFT
	  fftw_execute(ws.p_forw_src);

	  // Compute the element-wise product on the packed terms
	  // Let's put the element wise products in ws.in_kernel
	  double re_s, im_s, re_k, im_k;
	  for(ptr = ws.out_src, ptr2 = ws.out_kernel, ptr3 = ws.src_fft, ptr_end = ws.out_src+2*ws.h_fftw * (ws.w_fftw/2+1); ptr != ptr_end ; ptr+=2, ptr2+=2, ptr3+=2)
	    {
	      re_s = *ptr;
	      im_s = *(ptr+1);
	      re_k = *ptr2;
	      im_k = *(ptr2+1);
	      *ptr3 = re_s * re_k - im_s * im_k;
	      *(ptr3+1) = re_s * im_k + im_s * re_k;
	    }

	  // Compute the backward FFT
	  // Carefull, The backward FFT does not preserve the output
	  fftw_execute(ws.p_back);
	  // Scale the transform
	  unsigned int wh = ws.w_fftw*ws.h_fftw;
	  double wh_ = wh;
	  for(ptr = ws.dst_fft, ptr_end = ws.dst_fft + wh ; ptr != ptr_end ; ++ptr)
	    *ptr /= wh_;
	}

	void factorize (unsigned int n,
			unsigned int *n_factors,
			unsigned int factors[]) const {
	  unsigned int nf = 0;
	  unsigned int ntest = n;
	  unsigned int factor;
	  unsigned int i = 0;

	  if (n == 0)
	    throw std::invalid_argument("Length n must be positive integer");

	  if (n == 1) {
	    factors[0] = 1;
	    *n_factors = 1;
	    return ;
	  }

	  /* deal with the implemented factors */

	  while (FFTW_FACTORS[i] && ntest != 1)
	    {
	      factor = FFTW_FACTORS[i];
	      while ((ntest % factor) == 0)
		{
		  ntest = ntest / factor;
		  factors[nf] = factor;
		  nf++;
		}
	      i++;
	    }

	  // Ok that's it
	  if(ntest != 1) {
	    factors[nf] = ntest;
	    nf++;
	  }

	  /* check that the factorization is correct */
	  {
	    unsigned int product = 1;

	    for (i = 0; i < nf; i++)
	      product *= factors[i];

	    if (product != n)
	      throw std::runtime_error("factorization failed");
	  }

	  *n_factors = nf;
	}


	bool is_optimal(unsigned int n)  const {
	  // We check that n is not a multiple of 4*4*4*2
	  if(n % 4*4*4*2 == 0)
	    return false;

	  unsigned int nf=0;
	  unsigned int factors[64];
	  unsigned int i = 0;
	  factorize(n, &nf, factors);

	  // We just have to check if the last factor belongs to GSL_FACTORS
	  while(FFTW_FACTORS[i]) {
	    if(factors[nf-1] == FFTW_FACTORS[i])
	      return true;
	    ++i;
	  }
	  return false;
	}

	unsigned int find_closest_factor(unsigned int n) const {
	  unsigned int j = n;
	  while(!is_optimal(j)) ++j;
	  return j;
	}

	void init_workspace() {
	  ws.h_fftw = find_closest_factor(height + kernel_height/2);
	  ws.w_fftw = find_closest_factor(width  + kernel_width/2);

	  ws.in_src = new double[ws.h_fftw * ws.w_fftw];
	  ws.out_src = (double*) fftw_malloc(sizeof(fftw_complex) * ws.h_fftw * (ws.w_fftw/2+1));
	  ws.in_kernel = new double[ws.h_fftw * ws.w_fftw];
	  ws.out_kernel = (double*) fftw_malloc(sizeof(fftw_complex) * ws.h_fftw * (ws.w_fftw/2+1));

	  ws.src_fft = (double*) fftw_malloc(sizeof(fftw_complex) * ws.h_fftw * (ws.w_fftw/2+1));
	  ws.dst_fft = new double[ws.h_fftw * ws.w_fftw];
	  ws.dst = new double[height * width];

	  // Initialization of the plans
	  ws.p_forw_src = fftw_plan_dft_r2c_2d(ws.h_fftw, ws.w_fftw, ws.in_src, (fftw_complex*)ws.out_src, FFTW_ESTIMATE);
	  ws.p_forw_kernel = fftw_plan_dft_r2c_2d(ws.h_fftw, ws.w_fftw, ws.in_kernel, (fftw_complex*)ws.out_kernel, FFTW_ESTIMATE);

	  // The backward FFT takes ws.out_kernel as input !!
	  ws.p_back = fftw_plan_dft_c2r_2d(ws.h_fftw, ws.w_fftw, (fftw_complex*)ws.src_fft, ws.dst_fft, FFTW_ESTIMATE);
	}

	void clear_workspace()
	{
	  if(ws.in_src) delete[] ws.in_src;
	  ws.in_src = 0;
	  if(ws.out_src) fftw_free((fftw_complex*)ws.out_src);
	  ws.out_src = 0;
	  if(ws.in_kernel) delete[] ws.in_kernel;
	  ws.in_kernel = 0;
	  if(ws.out_kernel) fftw_free((fftw_complex*)ws.out_kernel);
	  ws.out_kernel = 0;

	  if(ws.src_fft) fftw_free((fftw_complex*)ws.src_fft);
	  ws.src_fft = 0;
	  if(ws.dst_fft) delete[] ws.dst_fft;
	  ws.dst_fft = 0;
	  if(ws.dst) delete[] ws.dst;
	  ws.dst = 0;

	  // Destroy the plans
	  if(ws.p_forw_src) fftw_destroy_plan(ws.p_forw_src);
	  ws.p_forw_src = 0;
	  if(ws.p_forw_kernel) fftw_destroy_plan(ws.p_forw_kernel);
	  ws.p_forw_kernel = 0;
	  if(ws.p_back) fftw_destroy_plan(ws.p_back);
	  ws.p_back = 0;
	}

      public:
	void convolve() const {

	  // Compute the circular convolution
	  fftw_circular_convolution();

	  // Depending on the type of convolution one is looking for, we extract the appropriate part of the result from out_src
	  unsigned int h_offset, w_offset;

	  // Same linear convolution
	  // Here we just keep the first [h_filt/2:h_filt/2+h_dst-1 ; w_filt/2:w_filt/2+w_dst-1] real part elements of out_src
	  h_offset = kernel_height/2;
	  w_offset = kernel_width/2;
	  for(unsigned int i = 0 ; i < height ; ++i) {
	    double* dest      = &ws.dst[i*width];
	    double* src_begin = &ws.dst_fft[(i+h_offset)*ws.w_fftw+w_offset];
	    double* src_end   = src_begin + width;
	    std::copy(src_begin,src_end,dest);
	  }
	}

      public:

	enum class KernelType : char{
	  Gaussian,
	    Linear
	    };

   
	Convolution(unsigned int w, unsigned int h, double sigma, KernelType kernel_type) 
	  : width(w), height(h), size(w*h), data(w*h),
	    kernel_width(2*width -1), kernel_height(2*height -1),
	    kernel(new double[kernel_width*kernel_height]),
	    FFTW_FACTORS({{13,11,7,5,3,2,0}}) {
	  init_workspace();
	  double mid_w = kernel_width/2;
	  double mid_h = kernel_height/2;
	  double* k;
	  double hh,ww;
	  double norm = 0;

	  double sigma2, eh, ew;
	  switch(kernel_type) {
	  case KernelType::Gaussian:
	    sigma2 = sigma*sigma;
	    for(hh = 0, k = kernel; hh < kernel_height; ++hh) {
              eh = exp(-(hh-mid_h)*(hh-mid_h)/(2.0*sigma2));
              for(ww = 0; ww < kernel_width; ++ww, ++k) {
		ew = exp(-(ww-mid_w)*(ww-mid_w)/(2.0*sigma2));
		*k = eh*ew;
		norm += *k;
              }
	    }
	    break;
	  case KernelType::Linear:
	    for(hh = 0, k = kernel; hh < kernel_height; ++hh) {
	      if(fabs(hh - mid_h) > sigma)
                for(ww = 0; ww < kernel_width; ++ww, ++k)
		  *k = 0;
	      else {
                eh = 1.0 - fabs(hh - mid_h)/sigma;
                for(ww = 0; ww < kernel_width; ++ww, ++k) {
		  if(fabs(ww - mid_w) > sigma)
		    *k = 0;
		  else {
		    ew = 1.0 - fabs(ww - mid_w)/sigma; 
		    *k = eh * ew;
		    norm += *k;
		  }

                }
	      } 
	    }
	    break;
	  default:
	    throw std::invalid_argument("Bad KernelType");
	  }

	  double* k_end = kernel + (kernel_height*kernel_width);
	  for(k = kernel; k < k_end; ++k) *k /= norm;
	  fftw_compute_kernel();
	}

	virtual ~Convolution() {
	  delete [] kernel;
	  clear_workspace();
	}

	void clear() {
	  std::fill(data.begin(), data.end(), 0);
	}
      };
    }
  }


  namespace tab1d {
    namespace fft {
      class Convolution : public xsom::tab::fft::Convolution {
      public:

	const Mapping& mapping;

      public:

	Convolution(double sigma, KernelType kernel_type, const Mapping& m) :
	  xsom::tab::fft::Convolution(m.size, 1, sigma, kernel_type),
	  mapping(m) {}

	virtual ~Convolution() {}

	void learn(std::function<double (double)> fct_value_at) {

	  unsigned int idx;
	  auto c   = this->data.begin();
	  for(idx = 0; idx < this->width; ++idx, ++c) {
	    double pos = mapping.index2pos(idx);
	    *c = fct_value_at(pos);
	  }
	}

	void fill_plot(std::vector<ccmpl::Point>& points) const {
	  points.clear();

	  unsigned int idx;
	  double* c = this->ws.dst;
	  for(idx = 0; idx < this->width; ++idx, ++c) {
	    double pos = mapping.index2pos(idx);
	    points.push_back({pos,*c});
	  }
	}
	
	void fill_plot_noconv(std::vector<ccmpl::Point>& points) const {
	  points.clear();

	  unsigned int idx;
	  auto c   = this->data.begin();
	  for(idx = 0; idx < this->width; ++idx, ++c) {
	    double pos = mapping.index2pos(idx);
	    points.push_back({pos,*c});
	  }
	}

	double operator()(double pos) const {
	  return this->ws.dst[mapping.pos2rank(pos)];
	}
	
	double get_noconv(double pos) const {
	  return this->data[mapping.pos2rank(pos)];
	}
	
	double bmu() const {
	  this->convolve();
	  unsigned int w;
	  double* it = this->ws.dst;
	  unsigned int wmax = 0;
	  double       max  = *it;
	  double       m;
	  for(w=0;w<this->width;++w,++it) 
	    if((m = *it) > max) {
	      wmax = w;
	      max = m;
	    }
	  return mapping.index2pos(wmax);
	}
	//    double bmu_stochastic() const {
	//        this->convolve();
	//        unsigned int w;
	//        double* it = this->ws.dst;
	//        double       sum  = 0;
	//        for(w=0; w < this->width; ++w,++it)
	//            sum += *it;
	//
	//        double p = std::rand()/double(RAND_MAX) * sum;
	//        it = this->ws.dst;
	//        sum = *it;
	//        w = 0;
	//        while(sum < p) {
	//            ++w;
	//            sum += *(it++);
	//        }
	//
	//        return mapping.index2pos(w);
	//    }
      };

      inline Convolution convolution(double sigma, Convolution::KernelType kernel_type,
				     const Mapping& m) {
	return Convolution(sigma, kernel_type, m);
      }
    }
  }



  namespace tab2d {
    namespace fft {
      class Convolution : public xsom::tab::fft::Convolution {
      private:

	std::function<bool (const xsom::Point2D<double>&)> pos_is_valid;
      
      public:
	const Mapping& mapping;
      
	template<typename fctPOS_IS_VALID>
	Convolution(double sigma, KernelType kernel_type,
		    const Mapping& m,
		    const fctPOS_IS_VALID& fct_pos_is_valid) 
	  : xsom::tab::fft::Convolution(m.size.w,m.size.h,sigma,kernel_type),
	    pos_is_valid(fct_pos_is_valid),
	    mapping(m) {
	}

	virtual ~Convolution() {}     
		     
	virtual void learn(std::function<double (const xsom::Point2D<double>&)> fct_value_at) {
	  
	  xsom::Index2D idx;
	  auto c   = this->data.begin();
	  for(idx.h = 0; idx.h < this->height; ++idx.h)
	    for(idx.w = 0; idx.w < this->width; ++idx.w, ++c) {
	      xsom::Point2D<double> pos = mapping.index2pos(idx);
	      if(pos_is_valid(pos))
		*c = fct_value_at(pos);
	      else
		*c = 0;
	    }
	}
      
	virtual void fill_plot(std::vector<ccmpl::ValueAt>& points) const {
	  points.clear();

	  xsom::Index2D idx;
	  double* c = this->ws.dst;
	  for(idx.h = 0; idx.h < this->height; ++idx.h)
	    for(idx.w = 0; idx.w < this->width; ++idx.w, ++c) {
	      xsom::Point2D<double> pos = mapping.index2pos(idx);
	      if(pos_is_valid(pos))
		points.push_back({pos.x,pos.y,*c});
	    }
	}
      
      
	virtual void fill_plot_noconv(std::vector<ccmpl::ValueAt>& points) const {
	  points.clear();

	  xsom::Index2D idx;
	  auto c   = this->data.begin();
	  for(idx.h = 0; idx.h < this->height; ++idx.h)
	    for(idx.w = 0; idx.w < this->width; ++idx.w, ++c) {
	      xsom::Point2D<double> pos = mapping.index2pos(idx);
	      if(pos_is_valid(pos))
		points.push_back({pos.x,pos.y,*c});
	    }
	}

	double operator()(const xsom::Point2D<double>& pos) const {
	  return this->ws.dst[mapping.pos2rank(pos)];
	}
	
	double get_noconv(const xsom::Point2D<double>& pos) const {
	  return this->data[mapping.pos2rank(pos)];
	}
      
	virtual xsom::Point2D<double> bmu() const {
	  this->convolve();
	  unsigned int w,h;
	  double* it = this->ws.dst;
	  unsigned int wmax = 0;
	  unsigned int hmax = 0;
	  double       max  = *it;
	  double       m;
	  for(h=0;h<this->height;++h)
	    for(w=0;w<this->width;++w,++it) 
	      if((m = *it) > max) {
		wmax = w;
		hmax = h;
		max = m;
	      }
	  return mapping.index2pos({wmax,hmax});
	}
      };

      template<typename fctPOS_IS_VALID> 
      Convolution convolution(double sigma, Convolution::KernelType kernel_type,
			      const Mapping& m,
			      const fctPOS_IS_VALID& fct_pos_is_valid) {
	return Convolution(sigma,kernel_type, m,fct_pos_is_valid);
      }
    }
  }

  

  
}
