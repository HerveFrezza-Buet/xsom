#pragma once

#include <fftw3.h>
#include <vector>
#include <array>
#include <cmath>
#include <random>

#include <xsomTable.hpp>
#include <xsomUtility.hpp>

namespace xsom {
  
  namespace tab {
    namespace fft {
      enum class KernelType : char{
	Gaussian,
	  Linear
	  };

	  enum class PaddingType: char {
		  Zero,
		  Constant
	  };

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
	std::vector<double>& data;
	
      public:
	mutable Workspace ws;
      protected:
	unsigned int kernel_width,kernel_height; 
	double* kernel;
	std::array<unsigned int,7> FFTW_FACTORS;
	PaddingType padding_type;

	void fftw_compute_kernel() {
	  double * ptr, *ptr_end;
	  for(ptr = ws.in_kernel, ptr_end = ws.in_kernel + ws.h_fftw*ws.w_fftw ; ptr != ptr_end ; ++ptr)
	    *ptr = 0.0;

	  for(unsigned int i = 0 ; i < kernel_height ; ++i)
	    for(unsigned int j = 0 ; j < kernel_width ; ++j)
	      ws.in_kernel[i*ws.w_fftw+j] = kernel[i*kernel_width + j];

	  fftw_execute(ws.p_forw_kernel);
	}

	void fftw_circular_convolution() const {
	  double * ptr, *ptr_end, *ptr2, *ptr3;

	  // Reset the content of ws.in_src
	  std::fill(ws.in_src, ws.in_src + ws.h_fftw*ws.w_fftw, 0.0);

	  // Then we build our periodic signals
	  auto in_src_line     = ws.in_src;
	  auto in_src_end      = ws.in_src + height*ws.w_fftw;
	  auto data_line       = std::data(data);
	  auto data_line_end   = data_line + width;
	  for(; in_src_line != in_src_end; in_src_line += ws.w_fftw, data_line = data_line_end, data_line_end += width)
	    std::copy(data_line, data_line_end, in_src_line);

	  // This is slower that the previous lines.
	  // for(unsigned int i = 0 ; i < height ; ++i)
	  //   for(unsigned int j = 0 ; j < width ; ++j)
	  //     ws.in_src[i*ws.w_fftw+j] = data[i*width + j];
	  
	  if(this->padding_type == PaddingType::Constant) {
		  // In this case, we need to duplicate the border values of the
		  // signal in the right regions of the extended w_fftw
		  // ----------------------
		  // |tl t1 t2  .... tn tr|
		  // |l1                r1|
		  // |l2                r2|
		  // |..                ..|
		  // |lm                rm|
		  // |bl b1 b2  .... bn br|
		  // ----------------------
		  //
		  // kh = int(kernel_height/2)
		  // kw = int(kernel_width/2)
		  //
		  //  ---------------------- 1  2 3  .. kw         1  2  3  .. kw
		  //  |tl t1 t2  .... tn tr|tr tr tr .. tr  0 .. 0 tl tl tl .. tl
		  //  |l1                r1|r1 r1 r1 .. r1  0 .. 0 l1 l1 l1 .. l1
		  //  |l2                r2|r2 r2 r2 .. r2  0 .. 0 l2 l2 l2 .. l2
		  //  |..                ..|
		  //  |lm                rm|rm rm rm .. rm  0 .. 0 lm lm lm .. lm
		  //  |bl b1 b2  .... bn br|br br br .. br  0 .. 0 bl bl bl .. bl
		  //  ----------------------
		  //1  bl b1 b2  .... bn br br br br .. br  0 .. 0 bl bl bl .. bl
		  //2  bl b1 b2  .... bn br br br br .. br  0 .. 0 bl bl bl .. bl 
		  //3  bl b1 b2  .... bn br br br br .. br  0 .. 0 bl bl bl .. bl
		  //.  .. .. ..  .... .. .. .. .. .. .. ..  0 .. 0 .. .. .. .. ..
		  //kh bl b1 b2  .... bn br br br br .. br  0 .. 0 bl bl bl .. bl
		  //
		  //   0  0  0   .... 0  0  0  0  0  .. 0   0 .. 0 0  0  0  .. 0 
		  //   .. .. ..  .... .. .. .. .. .. .. ..  0 .. 0 .. .. .. .. ..
		  //   0  0  0   .... 0  0  0  0  0  .. 0   0 .. 0 0  0  0  .. 0 
		  //1  tl t1 t2  .... tn tr tr tr tr .. tr  0 .. 0 tl tl tl .. tl
		  //2  tl t1 t2  .... tn tr tr tr tr .. tr  0 .. 0 tl tl tl .. tl 
		  //3  tl t1 t2  .... tn tr tr tr tr .. tr  0 .. 0 tl tl tl .. tl 
		  //.  .. .. ..  .... .. .. .. .. .. .. ..  0 .. 0 .. .. .. .. ..
		  //kh tl t1 t2  .... tn tr tr tr tr .. tr  0 .. 0 tl tl tl tl tl

		  int kh = int(kernel_height/2);
		  int kw = int(kernel_width/2);
		  double* first_element_ptr;

		  // Duplicate the rightmost colum
		  for(int i = 0 ; i <= kh; ++i)
			  std::fill(ws.in_src + i * ws.w_fftw + width, 
					    ws.in_src + i * ws.w_fftw + width + kw, 
						data[i*width + width - 1]);
		  
		  // Copy the first column on the top right of in_src
		  for(int i = 0; i <= kh; ++i)
			  std::fill(ws.in_src + i * ws.w_fftw + ws.w_fftw - kw, 
					    ws.in_src + i * ws.w_fftw + ws.w_fftw, 
						data[i * width]);

		  // Duplicate the bottom line
		  for(int i = 0; i <= kh; ++i)
			  std::copy(ws.in_src + (height-1)*ws.w_fftw, 
					    ws.in_src + (height-1)*ws.w_fftw + width,
						ws.in_src + (height+i)*ws.w_fftw);

		  // Make the block full of  'br'
		  first_element_ptr = ws.in_src + height * ws.w_fftw + width;
		  for(int i = 0; i <= kh; ++i)
			  std::fill(first_element_ptr + i * ws.w_fftw, 
					    first_element_ptr + i * ws.w_fftw + kw, 
						data[(height-1)*width + width-1]);
		 
		  // make the block of bl
		  first_element_ptr = ws.in_src + height * ws.w_fftw + ws.w_fftw - 1 - kw;
		  for(int i = 0; i <= kh; ++i)
			  std::fill(first_element_ptr + i * ws.w_fftw, 
					    first_element_ptr + i * ws.w_fftw + kw, 
						data[(height-1)*width]);
		  
		  // Duplicate the top line
		  first_element_ptr = &(data[0]);
		  for(int i = 0; i < kh; ++i)
			  std::copy(first_element_ptr,
					  first_element_ptr + width, 
					  ws.in_src + (ws.h_fftw - 1) * ws.w_fftw 
					            - (kh - 1) * ws.w_fftw + i * ws.w_fftw);
 
		  // Make the block of 'tr'
		  first_element_ptr = ws.in_src + (ws.h_fftw - 1) * ws.w_fftw
			                 - (kh - 1) * ws.w_fftw + width;
		  for(int i = 0; i < kh; ++i)
			  std::fill(first_element_ptr + i * ws.w_fftw, 
					    first_element_ptr + i * ws.w_fftw + kw, 
						data[width-1]);

		  // Make the block full of  'tl'
		  first_element_ptr = ws.in_src + (ws.h_fftw - 1) * ws.w_fftw
			                 - (kh - 1) * ws.w_fftw + ws.w_fftw - kw;
		  for(int i = 0; i < kh; ++i)
			  std::fill(first_element_ptr + i * ws.w_fftw, 
					    first_element_ptr + i * ws.w_fftw + kw, 
						data[0]);
	  }


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
	  double inv_wh_ = 1.0/wh;
	  for(ptr = ws.dst_fft, ptr_end = ws.dst_fft + wh ; ptr != ptr_end ; ++ptr)
	    *ptr *= inv_wh_;
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
	  if(this->padding_type == PaddingType::Zero) {
		  // We can extend by half the kernel size because it's just
		  // zero padding so the left part of the kernel will overlap the 
		  // same zeros than the right part of the kernel
		  ws.h_fftw = find_closest_factor(height + kernel_height/2);
		  ws.w_fftw = find_closest_factor(width  + kernel_width/2);
	  }
	  else if(this->padding_type == PaddingType::Constant) {
		  // For this case, we may need to extend by kernel_height
		  // and kernel_width, not halved. 
		  ws.h_fftw = find_closest_factor(height + kernel_height);
		  ws.w_fftw = find_closest_factor(width  + kernel_width);
	  }
	  
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

   
	Convolution(unsigned int w, unsigned int h, double sigma,
		    std::vector<double>& data,
		    KernelType kernel_type,
			PaddingType padding_type) 
	  : width(w), height(h), size(w*h),
	    data(data),
	    kernel_width(2*width -1), kernel_height(2*height -1),
	    kernel(new double[kernel_width*kernel_height]),
	    FFTW_FACTORS({{13,11,7,5,3,2,0}}),
	   padding_type(padding_type)	{
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
              eh = std::exp(-(hh-mid_h)*(hh-mid_h)/(2.0*sigma2));
              for(ww = 0; ww < kernel_width; ++ww, ++k) {
		ew = std::exp(-(ww-mid_w)*(ww-mid_w)/(2.0*sigma2));
		*k = eh*ew;
		norm += *k;
              }
	    }
	    break;
	  case KernelType::Linear:
	    for(hh = 0, k = kernel; hh < kernel_height; ++hh) {
	      if(std::fabs(hh - mid_h) > sigma)
                for(ww = 0; ww < kernel_width; ++ww, ++k)
		  *k = 0;
	      else {
                eh = 1.0 - std::fabs(hh - mid_h)/sigma;
                for(ww = 0; ww < kernel_width; ++ww, ++k) {
		  if(std::fabs(ww - mid_w) > sigma)
		    *k = 0;
		  else {
		    ew = 1.0 - std::fabs(ww - mid_w)/sigma; 
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
      class Convolution : public xsom::tab1d::Table<double> {
      public :

	mutable xsom::tab::fft::Convolution convolution;

      public:

	Convolution(const xsom::tab1d::Mapping& m, double sigma, xsom::tab::fft::KernelType kernel_type, xsom::tab::fft::PaddingType padding_type)
	  : xsom::tab1d::Table<double>(m),
	  convolution(m.size, 1, sigma, this->content, kernel_type, padding_type) {}
	
	template<typename fctPOS_IS_VALID>
	Convolution(const xsom::tab1d::Mapping& m, double sigma, xsom::tab::fft::KernelType kernel_type, xsom::tab::fft::PaddingType padding_type,
		    const fctPOS_IS_VALID& fct_pos_is_valid)
	  : xsom::tab1d::Table<double>(m, fct_pos_is_valid),
	  convolution(m.size, 1, sigma, this->content, kernel_type, padding_type) {}
	
	double get(double pos) const {
	  return convolution.ws.dst[this->mapping.pos2rank(pos)];
	}
      
	double get_noconv(double pos) const {
	  return this->xsom::tab1d::Table<double>::get(pos);
	}


	double operator()(double pos) const {
	  return get(pos);
	}
	
	void convolve() {
	  convolution.convolve();
	}	

	void learn(std::function<double (double)> fct_value_at) {
	  this->xsom::tab1d::Table<double>::learn(fct_value_at);
	}
	
	double bmu() const {
	  convolution.convolve();
	  
	  if(this->pos_is_valid) {
	    unsigned int length  = this->mapping.length;
	    auto c               = convolution.ws.dst;
	    double best_pos      = 0;
	    double best_value    = 0;
	    unsigned int rank    = 0;
	  
	    // Let us find the first valid data.
	    for(rank = 0; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos)) {
		best_pos = pos;
		best_value = *c;
		break;
	      }
	    }

	    for(++rank,++c; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos) && (best_value < *c)) {
		best_pos = pos;
		best_value = *c;
	      }
	    }

	    return best_pos;
	  }
	  else {
	    auto begin = convolution.ws.dst;
	    auto end   = begin + this->mapping.length;
	    return this->mapping.rank2pos((unsigned int)(std::distance(begin, std::max_element(begin, end))));
	  }
	}

	
	template<typename RANDOM_DEVICE>
	double random_bmu(RANDOM_DEVICE& rd) const {
	  convolution.convolve();
	  
	  unsigned int length            = this->mapping.length;
	  auto c                         = convolution.ws.dst;
	  std::vector<double> best_pos;
	  double              best_value = 0;
	  unsigned int        rank       = 0;
	
	  if(pos_is_valid) {
	  
	    // Let us find the first valid data.
	    for(rank = 0; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos)) {
		best_pos.push_back(pos);
		best_value = *c;
		break;
	      }
	    }

	    for(++rank,++c; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos)) {
		if(best_value < *c) {
		  best_pos.clear();
		  best_pos.push_back(pos);
		  best_value = *c;
		}
		else if(best_value == *c)
		  best_pos.push_back(pos);
	      }
	    }

	  }
	  else {
	    // Let us find the first valid data... 0 here.
	    best_pos.push_back(mapping.rank2pos(0));
	    best_value = *(c++);
	  
	    for(rank = 1; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(best_value < *c) {
		best_pos.clear();
		best_pos.push_back(pos);
		best_value = *c;
	      }
	      else if(best_value == *c)
		best_pos.push_back(pos);
	    }
	  }

	  return best_pos[std::uniform_int_distribution<std::size_t>(0, best_pos.size()-1)(rd)];
	}
	
	
	double bmu_noconv() const {
	  return this->xsom::tab1d::Table<double>::bmu();
	}
	
	template<typename RANDOM_DEVICE>
	double random_bmu_noconv(RANDOM_DEVICE& rd) const {
	  return this->xsom::tab1d::Table<double>::random_bmu(rd);
	}

	void fill_line(std::vector<ccmpl::Point>& points) const {
	  points.clear();

	  unsigned int rank;
	  unsigned int length  = this->mapping.length;
	  auto c   = convolution.ws.dst;
	  for(rank = 0; rank < length; ++rank) {
	    auto pos = this->mapping.rank2pos(rank);
	    points.push_back({pos, *(c++)});
	  }
	}
	
	void fill_line_noconv(std::vector<ccmpl::Point>& points) const {
	  this->xsom::tab1d::Table<double>::fill_line(points);
	}


      };

      inline Convolution convolution(const xsom::tab1d::Mapping& m, double sigma, xsom::tab::fft::KernelType kernel_type, xsom::tab::fft::PaddingType padding_type) {
	return Convolution(m, sigma, kernel_type, padding_type);
      }

      template<typename fctPOS_IS_VALID>
      inline Convolution convolution(const xsom::tab1d::Mapping& m, double sigma, xsom::tab::fft::KernelType kernel_type, xsom::tab::fft::PaddingType padding_type, 
				     const fctPOS_IS_VALID& fct_pos_is_valid) {
	return Convolution(m, sigma, kernel_type, padding_type, fct_pos_is_valid);
      }
    }
  }



  namespace tab2d {
    namespace fft {

      class Convolution : public xsom::tab2d::Table<double> {
      public:

	mutable xsom::tab::fft::Convolution convolution;

      public:

	Convolution(const  xsom::tab2d::Mapping& m,
		    double sigma, xsom::tab::fft::KernelType kernel_type, 
			xsom::tab::fft::PaddingType padding_type)
	  : xsom::tab2d::Table<double>(m),
	  convolution(m.size.w, m.size.h, sigma, this->content, kernel_type, padding_type) {}
	
	template<typename fctPOS_IS_VALID>
	Convolution(const  xsom::tab2d::Mapping& m,
		    double sigma, xsom::tab::fft::KernelType kernel_type, xsom::tab::fft::PaddingType padding_type,
		    const fctPOS_IS_VALID& fct_pos_is_valid)
	  : xsom::tab2d::Table<double>(m,fct_pos_is_valid),
	  convolution(m.size.w, m.size.h, sigma, this->content, kernel_type, padding_type) {}


	void convolve() {
	  convolution.convolve();
	}
	
	double get_noconv(xsom::Point2D<double> pos) const {
	  return this->xsom::tab2d::Table<double>::get(pos);
	}
	
	double get(xsom::Point2D<double> pos) const {
	  return convolution.ws.dst[this->mapping.pos2rank(pos)];
	}
	
	double operator()(xsom::Point2D<double> pos) const {
	  return get(pos);
	}
	
	void learn(std::function<double (xsom::Point2D<double>)> fct_value_at) {
	  this->xsom::tab2d::Table<double>::learn(fct_value_at);
	}
	
	xsom::Point2D<double> bmu() const {
	  convolution.convolve();
	  
	  if(this->pos_is_valid) {
	    unsigned int          length    = this->mapping.length;
	    auto                  c         = convolution.ws.dst;
	    xsom::Point2D<double> best_pos;
	    double best_value               = 0;
	    unsigned int rank               = 0;
	  
	    // Let us find the first valid data.
	    for(rank = 0; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos)) {
		best_pos = pos;
		best_value = *c;
		break;
	      }
	    }

	    for(++rank,++c; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos) && (best_value < *c)) {
		best_pos = pos;
		best_value = *c;
	      }
	    }

	    return best_pos;
	  }
	  else {
	    auto begin = convolution.ws.dst;
	    auto end   = begin + this->mapping.length;
	    return this->mapping.rank2pos((unsigned int)(std::distance(begin, std::max_element(begin, end))));
	  }
	}


	template<typename RANDOM_DEVICE>
	xsom::Point2D<double>  random_bmu(RANDOM_DEVICE& rd) const {
	  convolution.convolve();
	  
	  unsigned int length                           = this->mapping.length;
	  auto c                                        = convolution.ws.dst;
	  std::vector<xsom::Point2D<double>> best_pos;
	  double       best_value                       = 0;
	  unsigned int rank                             = 0;
	
	  if(pos_is_valid) {
	  
	    // Let us find the first valid data.
	    for(rank = 0; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos)) {
		best_pos.push_back(pos);
		best_value = *c;
		break;
	      }
	    }

	    for(++rank,++c; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(pos_is_valid(pos)) {
		if(best_value < *c) {
		  best_pos.clear();
		  best_pos.push_back(pos);
		  best_value = *c;
		}
		else if(best_value == *c)
		  best_pos.push_back(pos);
	      }
	    }

	  }
	  else {
	    // Let us find the first valid data... 0 here.
	    best_pos.push_back(mapping.rank2pos(0));
	    best_value = *(c++);	    
	  
	    for(rank = 1; rank < length; ++rank, ++c) {
	      auto pos = mapping.rank2pos(rank);
	      if(best_value < *c) {
		best_pos.clear();
		best_pos.push_back(pos);
		best_value = *c;
	      }
	      else if(best_value == *c)
		best_pos.push_back(pos);
	    }
	  }

	  return best_pos[std::uniform_int_distribution<std::size_t>(0, best_pos.size()-1)(rd)];
	}
	
	
	xsom::Point2D<double> bmu_noconv() const {
	  return this->xsom::tab2d::Table<double>::bmu();
	}
	
	
	template<typename RANDOM_DEVICE>
	xsom::Point2D<double> random_bmu_noconv(RANDOM_DEVICE& rd) const {
	  return this->xsom::tab2d::Table<double>::random_bmu(rd);
	}

	
	template<typename RANDOM_DEVICE>
	void fill_surface(unsigned int nb_triangulation_points,
			  RANDOM_DEVICE& rd,
			  std::vector<ccmpl::ValueAt>& points) const {
	  if(nb_triangulation_points != this->delaunay.size()) {
	    this->delaunay.clear();
	    this->delaunay.reserve(nb_triangulation_points);
	    auto out = std::back_inserter(this->delaunay);
	    unsigned int width  = this->mapping.size.w;
	    unsigned int height = this->mapping.size.h;
	    unsigned int nb=0;
	    if(this->pos_is_valid) {
	      while(nb < nb_triangulation_points) {
		auto idx = xsom::index2d(width, height, rd);
		auto pos = this->mapping.index2pos(idx);
		if(this->pos_is_valid(pos)) {
		  *(out++) = {this->mapping.index2rank(idx),pos};
		  ++nb;
		}
	      }
	    }
	    else 
	      while(nb < nb_triangulation_points) {
		auto idx = xsom::index2d(width, height, rd);
		auto pos = this->mapping.index2pos(idx);
		*(out++) = {this->mapping.index2rank(idx),pos};
		++nb;
	      }
	      
	  }
	  
	  points.clear();
	  for(auto& ip : this->delaunay)
	    points.push_back({ip.second.x,ip.second.y,this->convolution.ws.dst[ip.first]});
	}
	

	template<typename RANDOM_DEVICE>
	void fill_surface_noconv(unsigned int nb_triangulation_points,
				 RANDOM_DEVICE& rd,
				 std::vector<ccmpl::ValueAt>& points) const {
	  this->xsom::tab2d::Table<double>::fill_surface(nb_triangulation_points,points,rd);
	}

	void fill_image_gray(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
			     unsigned int& width, unsigned int& depth) const {
	  depth = 1;
	  width = this->mapping.size.w;
	  x.clear();
	  y.clear();
	  z.clear();
	  auto outx = std::back_inserter(x);
	  auto outy = std::back_inserter(y);
	  auto outz = std::back_inserter(z);
	
	  xsom::Index2D idx;
	  unsigned int height = this->mapping.size.h;
	  auto c = convolution.ws.dst;

	  idx.h = 0;
	  *(outy++) = this->mapping.index2pos({0,0}).y;
	  for(idx.w = 0; idx.w < width; ++idx.w) {
	    auto pos = this->mapping.index2pos(idx);
	    *(outx++) = pos.x;
	    *(outz++) = *(c++);
	  }
	
	  for(++idx.h; idx.h < height; ++idx.h) {
	    idx.w = 0;
	    auto pos = this->mapping.index2pos(idx);
	    *(outy++) = pos.y;
	    *(outz++) = *(c++);
	    for(idx.w = 1; idx.w < width; ++idx.w)
	      *(outz++) = *(c++);
	  }
	}
	
	void fill_image_gray_noconv(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
				    unsigned int& width, unsigned int& depth) const {
	  return this->xsom::tab2d::Table<double>::fill_image_gray(x,y,z,width,depth);
	}
      };


      inline Convolution convolution(const xsom::tab2d::Mapping& m,
				     double sigma, xsom::tab::fft::KernelType kernel_type,
					 xsom::tab::fft::PaddingType padding_type) {
	return Convolution(m,sigma,kernel_type, padding_type);
      }
      
      template<typename fctPOS_IS_VALID> 
      Convolution convolution(const xsom::tab2d::Mapping& m,
			      double sigma, xsom::tab::fft::KernelType kernel_type,
			      xsom::tab::fft::PaddingType padding_type,
				  const fctPOS_IS_VALID& fct_pos_is_valid) {
	return Convolution(m,sigma,kernel_type, padding_type, fct_pos_is_valid);
      }
    }
  }
}
