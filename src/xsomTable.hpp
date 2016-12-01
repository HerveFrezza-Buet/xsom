#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>

#include <ccmpl.hpp>
#include <xsomPoint.hpp>

namespace xsom {

  namespace tab {

    template<typename INDEX, typename POSITION>
    class Mapping {
    public:
      
      POSITION min,max;
      INDEX size;

      typedef INDEX    index_type;
      typedef POSITION position_type;

      Mapping(POSITION pos_min, 
	      POSITION pos_max,
	      INDEX    size) 
	: min(pos_min), max(pos_max), size(size) {}

      INDEX        pos2index(POSITION      pos)  const {return INDEX();}
      POSITION     index2pos(INDEX         idx)  const {return POSITION();}
      unsigned int pos2rank (POSITION      pos)  const {return 0;}
      POSITION     rank2pos (unsigned int rank)  const {return POSITION();}
    };

    template<typename CONTENT, typename MAPPING>
    class Table {
    public:
      
      const MAPPING& mapping;
      std::vector<CONTENT> content;

      Table(const MAPPING& m) 
	: mapping(m), 
	  content(m.size)  {}

      CONTENT operator()(typename MAPPING::position_type pos) const {
	return content[mapping.pos2rank(pos)];
      }
      
      void clear(const CONTENT& val) {
	std::fill(content.begin(), content.end(), val);
      }

      void learn(std::function<CONTENT (typename MAPPING::position_type)> fct_value_at) {
	unsigned int rank;
	unsigned int size  = mapping.size;
	auto c   = content.begin();
	for(rank = 0; rank < size; ++rank, ++c) {
	  auto pos = mapping.rank2pos(rank);
	  *c = fct_value_at(pos);
	}
      }
      
      typename MAPPING::position_type bmu() const {
	auto it = std::max_element(content.begin(), content.end());
	return mapping.rank2pos((unsigned int)(std::distance(content.begin(), std::max_element(content.begin(), content.end()))));
      }
	  
      void write(std::ostream& os) const {
	for(auto& val : content) os << val << '\n';
      }

      void read(std::istream& is) {
	char sep;
	for(auto& val : content) {
	  is >> val;
	  is.get(sep);
	}
      }
      
    };
  }


  
  namespace tab1d {
    

    class Mapping : public xsom::tab::Mapping<unsigned int,double> {
    public:
      typedef unsigned int index_type;
      typedef double       position_type;
      double coefx,coefw;

      Mapping(position_type pos_min, position_type pos_max, unsigned int nb_pos) 
	: xsom::tab::Mapping<unsigned int,double>(pos_min,pos_max,nb_pos) {
	coefw = (nb_pos-1.)/(max-min);
	coefx = (max-min)/(nb_pos-1.);
      }
      
      unsigned int pos2index(double pos)  const {
	return (unsigned int)((pos-min)*coefw+.5);
      }
      
      double index2pos(unsigned int idx)  const {
	return idx*coefx+min;
      }
      
      unsigned int pos2rank (double pos)  const {
	return (unsigned int)((pos-min)*coefw+.5);
      }
      
      double rank2pos (unsigned int rank)  const {
	return rank*coefx+min;
      }
    };

    Mapping mapping(double pos_min, double pos_max, unsigned int nb_pos) {
      return Mapping(pos_min,pos_max,nb_pos);
    }


    template<typename CONTENT>
    class Table_ : public xsom::tab::Table<CONTENT,Mapping> {
    public:
      const Mapping& mapping;

      Table_(const Mapping& mapping) : xsom::tab::Table<CONTENT,Mapping>(mapping) {}

    };

    template<typename CONTENT>
    class Table : public Table_<CONTENT> {
    public:

      using Table_<CONTENT>::Table_;
      
      template<typename ValueOf>
      void fill_line(const ValueOf& value_of_content,
		     std::vector<ccmpl::Point>& points) {
	points.clear();

	unsigned int rank;
	unsigned int size  = this->mapping.size;
	auto c   = this->content.begin();
	for(rank = 0; rank < size; ++rank, ++c) {
	  auto pos = this->mapping.rank2pos(rank);
	  points.push_back({pos, value_of_content(*c)});
	}
      }
    };

    template<>
    class Table<double> : public Table_<double> {
    public:

      using Table_<double>::Table_;
      
      void fill_line(std::vector<ccmpl::Point>& points) {
	points.clear();

	unsigned int rank;
	unsigned int size  = this->mapping.size;
	auto c   = this->content.begin();
	for(rank = 0; rank < size; ++rank, ++c) {
	  auto pos = this->mapping.rank2pos(rank);
	  points.push_back({pos, *c});
	}
      }
    };
    

    template<typename CONTENT>
    Table<CONTENT> table(const Mapping& m) {
      return Table<CONTENT>(m);
    }
  }

    
  namespace tab2d {
    class Mapping : public xsom::tab::Mapping<xsom::Index2D, xsom::Point2D<double> > {
    public:
      typedef xsom::Point2D<double>               position_type;
      typedef xsom::Index2D                       index_type;
      double coefx,coefy,coefw,coefh;
      
      Mapping(position_type pos_min, position_type pos_max, index_type size)
	: xsom::tab::Mapping<xsom::Index2D, xsom::Point2D<double> >(pos_min, pos_max, size) {
	coefw = (size.w-1.)/(max.x-min.x);
	coefh = (size.h-1.)/(max.y-min.y);
	coefx = (max.x-min.x)/(size.w-1.);
	coefy = (max.y-min.y)/(size.h-1.);
      }
      
      index_type pos2index(position_type pos)  const {
	return { (unsigned int)((pos.x-this->min.x)*coefw+.5), (unsigned int)((pos.y-this->min.y)*coefh+.5)};
      }
      
      position_type index2pos(index_type idx)  const {
	return {idx.w*coefx+this->min.x, idx.h*coefy+this->min.y};
      }
      
      unsigned int pos2rank (position_type pos)  const {
	return (unsigned int)((pos.y-this->min.y)*coefh+.5)*size.w+(unsigned int)((pos.x-this->min.x)*coefw+.5);
      }
      
      position_type rank2pos (unsigned int rank)  const {
	return index2pos({rank % size.w, rank / size.w});
      }
    };
    
    Mapping mapping(const xsom::Point2D<double>& pos_min,
		    const xsom::Point2D<double>& pos_max,
		    const xsom::Index2D& size) {
      return Mapping(pos_min,pos_max,size);
    }


    template<typename CONTENT>
    class Table_ : public xsom::tab::Table<CONTENT,Mapping> {
    public:
      const Mapping& mapping;
      std::function<bool (const xsom::Point2D<double>&)> pos_is_valid;

      template<typename fctPOS_IS_VALID>
      Table_(const Mapping& m,
	     const fctPOS_IS_VALID& fct_pos_is_valid) 
	: xsom::tab::Table<CONTENT,Mapping>(mapping),
	pos_is_valid(fct_pos_is_valid) {}

      void learn(std::function<CONTENT (const xsom::Point2D<double>&)> fct_value_at) {
	unsigned int rank;
	unsigned int size  = mapping.size;
	auto c   = this->content.begin();
	for(rank = 0; rank < size; ++rank, ++c) {
	  auto pos = mapping.rank2pos(rank);
	  if(pos_is_valid(pos))
	    *c = fct_value_at(pos);
	}
      }

      
      template<typename ColorOf>
      void fill_palette(const ColorOf& ccmplcolor_of_content,
			std::vector<ccmpl::ColorAt>& points) {
	points.clear();

	xsom::Index2D idx;
	unsigned int width  = mapping.size.w;
	unsigned int height = mapping.size.h;
	auto c   = this->content.begin();
	for(idx.h = 0; idx.h < height; ++idx.h)
	  for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	    xsom::Point2D<double> pos = mapping.index2pos(idx);
	    if(this->pos_is_valid(pos))
	      points.push_back({pos.x,pos.y,ccmplcolor_of_content(*c)});
	  }
      }
      
      template<typename ColorOf>
      void fill_image_rgb(const ColorOf& ccmplcolor_of_content,
			   std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
			  unsigned int& width, unsigned int& depth) {
	depth = 3;
	width = mapping.size.w;
	x.clear();
	y.clear();
	z.clear();
	auto outx = std::back_inserter(x);
	auto outy = std::back_inserter(y);
	auto outz = std::back_inserter(z);
	
	xsom::Index2D idx;
	unsigned int height = mapping.size.h;
	auto c   = this->content.begin();

	idx.h = 0;
	*(outy++) = mapping.index2pos({0,0}).y;
	for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	  auto pos = mapping.index2pos(idx);
	  *(outx++) = pos.x;
	  auto color = ccmplcolor_of_content(*c);
	  *(outz++) = c.r;
	  *(outz++) = c.g;
	  *(outz++) = c.b;
	}
	
	for(++idx.h; idx.h < height; ++idx.h) {
	  idx.w = 0;
	  auto pos = mapping.index2pos(idx);
	  *(outy++) = pos.y;
	  *(outz++) = value_of_content(*c);
	  for(idx.w = 1; idx.w < width; ++idx.w, ++c) {
	    auto color = ccmplcolor_of_content(*c);
	    *(outz++) = c.r;
	    *(outz++) = c.g;
	    *(outz++) = c.b;
	  }
	}
      }
      
    };

    template<typename CONTENT>
    class Table : public Table_<CONTENT> {
    public:

      using Table_<CONTENT>::Table_;
      
      template<typename ValueOf>
      void fill_surface(const ValueOf& value_of_content,
			std::vector<ccmpl::ValueAt>& points) {
	points.clear();

	xsom::Index2D idx;
	unsigned int width  = this->mapping.size.w;
	unsigned int height = this->mapping.size.h;
	auto c   = this->content.begin();
	for(idx.h = 0; idx.h < height; ++idx.h)
	  for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	    xsom::Point2D<double> pos = this->mapping.index2pos(idx);
	    if(this->pos_is_valid(pos))
	      points.push_back({pos.x,pos.y,value_of_content(*c)});
	  }
      }

      template<typename ValueOf>
      void fill_image_gray(const ValueOf& value_of_content,
			   std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
			   unsigned int& width, unsigned int& depth) {
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
	auto c   = this->content.begin();

	idx.h = 0;
	*(outy++) = this->mapping.index2pos({0,0}).y;
	for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	  auto pos = this->mapping.index2pos(idx);
	  *(outx++) = pos.x;
	  *(outz++) = value_of_content(*c);
	}
	
	for(++idx.h; idx.h < height; ++idx.h) {
	  idx.w = 0;
	  auto pos = this->mapping.index2pos(idx);
	  *(outy++) = pos.y;
	  *(outz++) = value_of_content(*c);
	  for(idx.w = 1; idx.w < width; ++idx.w, ++c)
	    *(outz++) = value_of_content(*c);
	}
      }
    };


    template<>
    class Table<double> : public Table_<double> {
    public:

      using Table_<double>::Table_;
      
      void fill_surface(std::vector<ccmpl::ValueAt>& points) {
	points.clear();

	xsom::Index2D idx;
	unsigned int width  = this->mapping.size.w;
	unsigned int height = this->mapping.size.h;
	auto c   = this->content.begin();
	for(idx.h = 0; idx.h < height; ++idx.h)
	  for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	    xsom::Point2D<double> pos = this->mapping.index2pos(idx);
	    if(this->pos_is_valid(pos))
	      points.push_back({pos.x,pos.y,*c});
	  }
      }

      void fill_image_gray(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
			   unsigned int& width, unsigned int& depth) {
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
	auto c   = this->content.begin();

	idx.h = 0;
	*(outy++) = this->mapping.index2pos({0,0}).y;
	for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	  auto pos = this->mapping.index2pos(idx);
	  *(outx++) = pos.x;
	  *(outz++) = *c;
	}
	
	for(++idx.h; idx.h < height; ++idx.h) {
	  idx.w = 0;
	  auto pos = this->mapping.index2pos(idx);
	  *(outy++) = pos.y;
	  *(outz++) = *c;
	  for(idx.w = 1; idx.w < width; ++idx.w, ++c)
	    *(outz++) = *c;
	}
      }
    };

    template<typename CONTENT,
	     typename fctPOS_IS_VALID>
    Table<CONTENT> table(const Mapping& m,
			 const fctPOS_IS_VALID& fct_pos_is_valid) {
      return Table<CONTENT>(m,fct_pos_is_valid);
    }

  }
    
    
}
