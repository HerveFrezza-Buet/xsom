#pragma once

#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>
#include <random>

#include <ccmpl.hpp>
#include <xsomPoint.hpp>
#include <xsomUtility.hpp>

namespace xsom {

  namespace tab {

    template<typename INDEX, typename POSITION>
    class Mapping {
    public:
      
      POSITION min,max;
      INDEX size;
      unsigned int length;

      typedef INDEX    index_type;
      typedef POSITION position_type;

      Mapping(POSITION pos_min, 
	      POSITION pos_max,
	      INDEX    size) 
	: min(pos_min), max(pos_max), size(size),
	  length(0) // this has to be set by subclasses.
      {}

      INDEX        pos2index  (POSITION      pos)  const {return    INDEX();}
      POSITION     index2pos  (INDEX         idx)  const {return POSITION();}
      unsigned int pos2rank   (POSITION      pos)  const {return          0;}
      POSITION     rank2pos   (unsigned int rank)  const {return POSITION();}
      unsigned int index2rank (INDEX         idx)  const {return          0;}
      INDEX        rank2index (unsigned int rank)  const {return    INDEX();}
    };

    template<typename CONTENT, typename MAPPING>
    class Table {
    public:
      
      const MAPPING& mapping;
      std::vector<CONTENT> content;
      std::function<bool (typename MAPPING::position_type pos)> pos_is_valid;

      friend std::ostream& operator<<(std::ostream& os, const Table& t) {
	for(auto d : t.content)
	  os << d << ' ';
	return os;
      }
      
      friend std::istream& operator>>(std::istream& is, Table& t) {
	for(auto& d : t.content)
	  is >> d;
	return is;
      }

      Table(const MAPPING& m) 
	: mapping(m), 
	  content(m.length),
	  pos_is_valid() {}
      
      template<typename fctPOS_IS_VALID>
      Table(const MAPPING& m,
	    const fctPOS_IS_VALID& fct_pos_is_valid) 
	: mapping(m), 
	  content(m.length),
	  pos_is_valid(fct_pos_is_valid) {}

      CONTENT get(typename MAPPING::position_type pos) const {
	return content[mapping.pos2rank(pos)];
      }
      
      CONTENT operator()(typename MAPPING::position_type pos) const {
	return get(pos);
      }
      
      void clear(const CONTENT& val) {
	std::fill(content.begin(), content.end(), val);
      }

      void learn(std::function<CONTENT (typename MAPPING::position_type)> fct_value_at) {
	unsigned int rank;
	unsigned int length  = mapping.length;
	auto c   = content.begin();

	if(pos_is_valid)
	  for(rank = 0; rank < length; ++rank, ++c) {
	    auto pos = this->mapping.rank2pos(rank);
	    if(pos_is_valid(pos))
	      *c = fct_value_at(pos);
	  }
	else
	  for(rank = 0; rank < length; ++rank) {
	    auto pos = mapping.rank2pos(rank);
	    *(c++) = fct_value_at(pos);
	  }
      }
      
      typename MAPPING::position_type bmu() const {
	if(pos_is_valid) {
	  unsigned int length  = mapping.length;
	  auto c               = content.begin();
	  typename MAPPING::position_type best_pos   =  typename MAPPING::position_type();
	  CONTENT                         best_value = CONTENT();
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
	else
	  return mapping.rank2pos((unsigned int)(std::distance(content.begin(), std::max_element(content.begin(), content.end()))));
      }

      /**
       * Thus function selects randomly among the best matching units (there may be several ones).
       */
      template<typename RANDOM_DEVICE>
      typename MAPPING::position_type random_bmu(RANDOM_DEVICE& rd) const {
	unsigned int length                                   = mapping.length;
	auto c                                                = content.begin();
	std::vector<typename MAPPING::position_type> best_pos;
	CONTENT                         best_value            = CONTENT();
	unsigned int rank                                     = 0;
	
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
	  
	  for(rank = 1;rank < length; ++rank, ++c) {
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

	return best_pos[std::uniform_int_distribution<std::size_t>(0,best_pos.size()-1)(rd)];
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
	length = size;
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
      
      unsigned int rank2index (unsigned int rank)  const {
	return rank;
      }
      
      unsigned int index2rank (unsigned int idx)  const {
	return idx;
      }
    };

    inline Mapping mapping(double pos_min, double pos_max, unsigned int nb_pos) {
      return Mapping(pos_min,pos_max,nb_pos);
    }


    template<typename CONTENT>
    class Table_ : public xsom::tab::Table<CONTENT,Mapping> {
    public:
      
      Table_(const Mapping& m)
	: xsom::tab::Table<CONTENT,Mapping>(m) {}
      
      template<typename fctPOS_IS_VALID>
      Table_(const Mapping& m,
	     const fctPOS_IS_VALID& fct_pos_is_valid)
	: xsom::tab::Table<CONTENT,Mapping>(m, fct_pos_is_valid) {}
    };

    template<typename CONTENT>
    class Table : public Table_<CONTENT> {
    public:

      Table(const Mapping& m)
	: Table_<CONTENT>(m) {}
      
      template<typename fctPOS_IS_VALID>
      Table(const Mapping& m,
	    const fctPOS_IS_VALID& fct_pos_is_valid)
	: Table_<CONTENT>(m, fct_pos_is_valid) {}

      void learn(std::function<CONTENT (double)> fct_value_at) {
	this->xsom::tab::Table<CONTENT,Mapping>::learn(fct_value_at);
      }
      
      double bmu() const {
	return this->xsom::tab::Table<CONTENT,Mapping>::bmu();
      }
      
      template<typename RANDOM_DEVICE>
      double random_bmu(RANDOM_DEVICE& rd) const {
	return this->xsom::tab::Table<CONTENT,Mapping>::random_bmu(rd);
      }
      
      template<typename ValueOf>
      void fill_line(const ValueOf& value_of_content,
		     std::vector<ccmpl::Point>& points) const {
	points.clear();

	unsigned int rank;
	unsigned int length  = this->mapping.length;
	auto c   = this->content.begin();
	for(rank = 0; rank < length; ++rank) {
	  auto pos = this->mapping.rank2pos(rank);
	  points.push_back({pos, value_of_content(*(c++))});
	}
      }
    };

    template<>
    class Table<double> : public Table_<double> {
    public:

      Table(const Mapping& m)
	: Table_<double>(m) {}
      
      template<typename fctPOS_IS_VALID>
      Table(const Mapping& m,
	    const fctPOS_IS_VALID& fct_pos_is_valid)
	: Table_<double>(m, fct_pos_is_valid) {}
      
      void learn(std::function<double (double)> fct_value_at) {
	this->xsom::tab::Table<double,Mapping>::learn(fct_value_at);
      }
      
      double bmu() const {
	return this->xsom::tab::Table<double,Mapping>::bmu();
      }
      
      template<typename RANDOM_DEVICE>
      double random_bmu(RANDOM_DEVICE& rd) const {
	return this->xsom::tab::Table<double,Mapping>::random_bmu(rd);
      }
      
      void fill_line(std::vector<ccmpl::Point>& points) const {
	points.clear();

	unsigned int rank;
	unsigned int length  = this->mapping.length;
	auto c   = this->content.begin();
	for(rank = 0; rank < length; ++rank) {
	  auto pos = this->mapping.rank2pos(rank);
	  points.push_back({pos, *(c++)});
	}
      }
    };
    

    template<typename CONTENT>
    Table<CONTENT> table(const Mapping& m) {
      return Table<CONTENT>(m);
    }
    
    template<typename CONTENT, typename fctPOS_IS_VALID>
    Table<CONTENT> table(const Mapping& m, const fctPOS_IS_VALID& fct_pos_is_valid) {
      return Table<CONTENT>(m,fct_pos_is_valid);
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
	length = size.w * size.h;
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

      unsigned int index2rank(index_type idx) const {
	return idx.w + size.w*idx.h;
      }
      
      index_type index2rank(unsigned int rank) const {
	return {rank % size.w, rank / size.w};
      }
    };
    
    inline Mapping mapping(const xsom::Point2D<double>& pos_min,
		    const xsom::Point2D<double>& pos_max,
		    const xsom::Index2D& size) {
      return Mapping(pos_min,pos_max,size);
    }


    template<typename CONTENT>
    class Table_ : public xsom::tab::Table<CONTENT,Mapping> {

    protected:

      mutable std::vector< std::pair<unsigned int,xsom::Point2D<double> > > delaunay;
      
    public:

      Table_(const Mapping& m) 
	: xsom::tab::Table<CONTENT,Mapping>(m),
	delaunay() {}
      
      template<typename fctPOS_IS_VALID>
      Table_(const Mapping& m,
	     const fctPOS_IS_VALID& fct_pos_is_valid) 
	: xsom::tab::Table<CONTENT,Mapping>(m,fct_pos_is_valid),
	delaunay() {}
      
      template<typename ColorOf, typename RANDOM_DEVICE>
      void fill_palette(const ColorOf& ccmplcolor_of_content,
			unsigned int nb_triangulation_points,
			RANDOM_DEVICE& rd,
			std::vector<ccmpl::ColorAt>& points) const {
	
	if(nb_triangulation_points != this->delaunay.size()) {
	  this->delaunay.clear();
	  this->delaunay.reserve(nb_triangulation_points);
	  auto out = std::back_inserter(this->delaunay);
	  unsigned int width  = this->mapping.size.w;
	  unsigned int height = this->mapping.size.h;
	  unsigned int nb=0;
	  if(this->pos_is_valid)
	    while(nb < nb_triangulation_points) {
	      auto idx = xsom::index2d(width, height, rd);
	      auto pos = this->mapping.index2pos(idx);
	      if(this->pos_is_valid(pos)) {
		*(out++) = {this->mapping.index2rank(idx),pos};
		++nb;
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
	for(auto& ip : this->delaunay) points.push_back({ip.second.x,ip.second.y,ccmplcolor_of_content(this->content[ip.first])});
      }
      
      
      template<typename ColorOf>
      void fill_image_rgb(const ColorOf& ccmplcolor_of_content,
			  std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
			  unsigned int& width, unsigned int& depth) const {
	depth = 3;
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
	for(idx.w = 0; idx.w < width; ++idx.w) {
	  auto pos = this->mapping.index2pos(idx);
	  *(outx++) = pos.x;
	  auto color = ccmplcolor_of_content(*(c++));
	  *(outz++) = color.r;
	  *(outz++) = color.g;
	  *(outz++) = color.b;
	}
	
	for(++idx.h; idx.h < height; ++idx.h) {
	  idx.w = 0;
	  auto pos = this->mapping.index2pos(idx);
	  *(outy++) = pos.y;
	  auto color = ccmplcolor_of_content(*(c++));
	  *(outz++) = color.r;
	  *(outz++) = color.g;
	  *(outz++) = color.b;
	  for(idx.w = 1; idx.w < width; ++idx.w) {
	    auto color = ccmplcolor_of_content(*(c++));
	    *(outz++) = color.r;
	    *(outz++) = color.g;
	    *(outz++) = color.b;
	  }
	}
      }
      
    };

    template<typename CONTENT>
    class Table : public Table_<CONTENT> {
    public:

      Table(const Mapping& m)
	: Table_<CONTENT>(m) {}
      
      template<typename fctPOS_IS_VALID>
      Table(const Mapping& m,
	    const fctPOS_IS_VALID& fct_pos_is_valid)
	: Table_<CONTENT>(m, fct_pos_is_valid) {}

      void learn(std::function<CONTENT (xsom::Point2D<double>)> fct_value_at) {
	this->xsom::tab::Table<CONTENT,Mapping>::learn(fct_value_at);
      }
      
      xsom::Point2D<double> bmu() const {
	return this->xsom::tab::Table<CONTENT,Mapping>::bmu();
      }
      
      template<typename RANDOM_DEVICE>
      xsom::Point2D<double> random_bmu(RANDOM_DEVICE& rd) const {
	return this->xsom::tab::Table<CONTENT,Mapping>::random_bmu(rd);
      }
      
      template<typename ColorOf, typename RANDOM_DEVICE>
      void fill_palette(const ColorOf& ccmplcolor_of_content,
			unsigned int nb_triangulation_points,
			RANDOM_DEVICE& rd,
			std::vector<ccmpl::ColorAt>& points) const {
	this->Table_<CONTENT>::fill_palette(ccmplcolor_of_content,
					    nb_triangulation_points,
					    rd,
					    points);
      }
      
      template<typename ColorOf>
      void fill_image_rgb(const ColorOf& ccmplcolor_of_content,
			  std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
			  unsigned int& width, unsigned int& depth) const {
	this->Table_<CONTENT>::fill_image_rgb(ccmplcolor_of_content,
					      x, y, z,
					      width, depth);
      }

      
      template<typename ValueOf, typename RANDOM_DEVICE>
      void fill_surface(const ValueOf& value_of_content,
			unsigned int nb_triangulation_points,
			RANDOM_DEVICE& rd,
			std::vector<ccmpl::ValueAt>& points) const {
	
	if(nb_triangulation_points != this->delaunay.size()) {
	  this->delaunay.clear();
	  this->delaunay.reserve(nb_triangulation_points);
	  auto out = std::back_inserter(this->delaunay);
	  unsigned int width  = this->mapping.size.w;
	  unsigned int height = this->mapping.size.h;
	  unsigned int nb=0;
	  if(this->pos_is_valid)
	    while(nb < nb_triangulation_points) {
	      auto idx = xsom::index2d(width, height, rd);
	      auto pos = this->mapping.index2pos(idx);
	      if(this->pos_is_valid(pos)) {
		*(out++) = {this->mapping.index2rank(idx),pos};
		++nb;
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
	for(auto& ip : this->delaunay) points.push_back({ip.second.x,ip.second.y,value_of_content(this->content[ip.first])});
      }

      template<typename ValueOf>
      void fill_image_gray(const ValueOf& value_of_content,
			   std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
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
	auto c   = this->content.begin();

	idx.h = 0;
	*(outy++) = this->mapping.index2pos({0,0}).y;
	for(idx.w = 0; idx.w < width; ++idx.w) {
	  auto pos = this->mapping.index2pos(idx);
	  *(outx++) = pos.x;
	  *(outz++) = value_of_content(*(c++));
	}
	
	for(++idx.h; idx.h < height; ++idx.h) {
	  idx.w = 0;
	  auto pos = this->mapping.index2pos(idx);
	  *(outy++) = pos.y;
	  *(outz++) = value_of_content(*(c++));
	  for(idx.w = 1; idx.w < width; ++idx.w)
	    *(outz++) = value_of_content(*(c++));
	}
      }
    };


    template<>
    class Table<double> : public Table_<double> {
    public:

      Table(const Mapping& m)
	: Table_<double>(m) {}
      
      template<typename fctPOS_IS_VALID>
      Table(const Mapping& m,
	    const fctPOS_IS_VALID& fct_pos_is_valid)
	: Table_<double>(m, fct_pos_is_valid) {}

      void learn(std::function<double (xsom::Point2D<double>)> fct_value_at) {
	this->xsom::tab::Table<double,Mapping>::learn(fct_value_at);
      }
      
      xsom::Point2D<double> bmu() const {
	return this->xsom::tab::Table<double,Mapping>::bmu();
      }
      
      template<typename RANDOM_DEVICE>
      xsom::Point2D<double> random_bmu(RANDOM_DEVICE& rd) const {
	return this->xsom::tab::Table<double,Mapping>::random_bmu(rd);
      }
      
      template<typename ColorOf, typename RANDOM_DEVICE>
      void fill_palette(const ColorOf& ccmplcolor_of_content,
			unsigned int nb_triangulation_points,
			RANDOM_DEVICE& rd,
			std::vector<ccmpl::ColorAt>& points) const {
	this->Table_<double>::fill_palette(ccmplcolor_of_content,
					   nb_triangulation_points,
					   rd,
					   points);
      }
      
      template<typename ColorOf>
      void fill_image_rgb(const ColorOf& ccmplcolor_of_content,
			  std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
			  unsigned int& width, unsigned int& depth) const {
	this->Table_<double>::fill_image_rgb(ccmplcolor_of_content,
					      x, y, z,
					      width, depth);
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
	  if(this->pos_is_valid)
	    while(nb < nb_triangulation_points) {
	      auto idx = xsom::index2d(width, height, rd);
	      auto pos = this->mapping.index2pos(idx);
	      if(this->pos_is_valid(pos)) {
		*(out++) = {this->mapping.index2rank(idx),pos};
		++nb;
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
	for(auto& ip : this->delaunay) points.push_back({ip.second.x,ip.second.y,this->content[ip.first]});
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
	auto c   = this->content.begin();

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
    };


    template<typename CONTENT>
    Table<CONTENT> table(const Mapping& m) {
      return Table<CONTENT>(m);
    }
    
    template<typename CONTENT, typename fctPOS_IS_VALID>
    Table<CONTENT> table(const Mapping& m, const fctPOS_IS_VALID& fct_pos_is_valid) {
      return Table<CONTENT>(m,fct_pos_is_valid);
    }

    constexpr unsigned int nb_lines(unsigned int nb_w, unsigned int nb_h) {
      return nb_w + nb_h + 2;
    }

    /**
     * This draw the 2D->2D mapping. It plots a grid, according to the
     * table coordinates. The position of the grid nodes depend on the
     * content of tha map.
     * @param NB_W, NB_H  Are the size of the displaied grid (i.e. nb of columns, nb of lines).
     * @param STEP Each line of the grid is built by dpanning one line of the table at every STEP elements.
     */
    template<unsigned int NB_W, unsigned int NB_H, unsigned int STEP>
    void fill_lines(const Table<xsom::Point2D<double> >& table, std::vector<std::vector<ccmpl::Point>>& lines) {
      lines.clear();
      lines.resize(nb_lines(NB_W,NB_H));
      auto it_line = lines.begin();
      auto begin   = table.content.begin();

      for(unsigned int w = 0; w <= NB_W; ++w, ++it_line) {
	it_line->clear();
	auto out = std::back_inserter(*it_line);
	unsigned int idx   = 0;
	auto         it    = begin;
	unsigned int ww    = 0;
	unsigned int bound = table.mapping.size.w-1;
	unsigned int idx0  = (unsigned int)(w*bound/(double)NB_W + .5);
	for(idx = idx0, it = begin + idx;
	    ww < bound;
	    idx +=  STEP*table.mapping.size.w, it = begin + idx, ww += STEP)
	  *(out++) = *it;
	*(out++) = *(begin + idx0 + bound*table.mapping.size.w);
      }

      for(unsigned int h = 0; h <= NB_H; ++h, ++it_line) {
	it_line->clear();
	auto out = std::back_inserter(*it_line);
	unsigned int idx   = 0;
	auto         it    = begin;
	unsigned int hh    = 0;
	unsigned int bound = table.mapping.size.h-1;
	unsigned int idx0  = (unsigned int)(h*bound/(double)NB_H + .5)*table.mapping.size.w;
	for(idx = idx0, it = begin + idx;
	    hh < bound;
	    idx += STEP, it = begin + idx, hh += STEP)
	  *(out++) = *it;
	*(out++) = *(begin + (idx0 + table.mapping.size.w-1));}
    }

    /**
     * Displays a line along the X axis of the 2D->2D projection. 
     * @param length_ratio 1 displays the whole x-axis. 
     * @param step The line is made of table elements taken at every step grid content.
     */
    inline void fill_x(const Table<xsom::Point2D<double> >& table, double length_ratio, unsigned int step, std::vector<ccmpl::Point>& points) {
      points.clear();
      auto         out   = std::back_inserter(points);
      auto         begin = table.content.begin();
      unsigned int bound = (unsigned int)(table.mapping.size.w*length_ratio+.5);
      unsigned int i     = 0;
      unsigned int idx   = 0;
      auto         it    = begin;
      for(; i < bound; idx += step, it = begin + idx, i += step)
	  *(out++) = *it;
    }


    /**
     * Displays a line along the Y axis of the 2D->2D projection. 
     * @param length_ratio 1 displays the whole y-axis. 
     * @param step The line is made of table elements taken at every step grid content.
     */
    inline void fill_y(const Table<xsom::Point2D<double> >& table, double length_ratio, unsigned int step, std::vector<ccmpl::Point>& points) {
      points.clear();
      auto         out   = std::back_inserter(points);
      auto         begin = table.content.begin();
      unsigned int bound = (unsigned int)(table.mapping.size.h*length_ratio+.5);
      unsigned int i     = 0;
      unsigned int idx   = 0;
      auto         it    = begin;
      for(; i < bound; idx += step*table.mapping.size.w, it = begin + idx, i += step)
	*(out++) = *it;
    }
  }
}
