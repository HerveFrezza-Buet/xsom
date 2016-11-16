#pragma once

#include <iostream>

namespace xsom {
  template<typename VALUE>
  struct Point2D {
    VALUE x,y;

    Point2D() : x(0), y(0) {}
    Point2D(const Point2D<VALUE>& p) : x(p.x), y(p.y) {}
    Point2D(VALUE xx,VALUE yy) : x(xx), y(yy) {}
    Point2D<VALUE>& operator=(const Point2D<VALUE>& p) {
      if(this != &p) {
	x = p.x;
	y = p.y;
      }
      return *this;
    }
    Point2D<VALUE> operator+(const Point2D<VALUE>& p) const {
      return Point2D<VALUE>(x+p.x,y+p.y);
    }
    Point2D<VALUE>& operator+=(const Point2D<VALUE>& p) {
      x+=p.x;
      y+=p.y;
      return *this;
    }
    Point2D<VALUE> operator-(const Point2D<VALUE>& p) const {
      return Point2D<VALUE>(x-p.x,y-p.y);
    }
    Point2D<VALUE>& operator-=(const Point2D<VALUE>& p) {
      x-=p.x;
      y-=p.y;
      return *this;
    }
    Point2D<VALUE> operator*(VALUE alpha) const {
      return Point2D<VALUE>(x*alpha,y*alpha);
    }
    Point2D<VALUE>& operator*=(VALUE alpha) {
      x*=alpha;
      y*=alpha;
      return *this;
    }
    Point2D<VALUE> operator/(VALUE alpha) const {
      return Point2D<VALUE>(x/alpha,y/alpha);
    }
    Point2D<VALUE>& operator/=(VALUE alpha) {
      x/=alpha;
      y/=alpha;
      return *this;
    }
    double operator*(const Point2D<VALUE>& p) const {
      return x*p.x+y*p.y;
    }
  };

  template<typename VALUE>
  std::ostream& operator<<(std::ostream& os,const Point2D<VALUE>& p) {
    os << p.x << ' ' << p.y;
    return os;
  }

  template<typename VALUE>
  std::istream& operator>>(std::istream& is, Point2D<VALUE>& p) {
    char sep;
    is >> p.x;
    is.get(sep);
    is >> p.y;
    return is;
  }
  
  
  template<typename VALUE>
  Point2D<VALUE> point2d(VALUE x,VALUE y) {
    return Point2D<VALUE>(x,y);
  }

  struct Index2D {
    unsigned int w,h;

    Index2D() : w(0), h(0) {}
    Index2D(const Index2D& p) : w(p.w), h(p.h) {}
    Index2D(unsigned int ww,unsigned int hh) : w(ww), h(hh) {}
    Index2D& operator=(const Index2D& p) {
      if(this != &p) {
	w = p.w;
	h = p.h;
      }
      return *this;
    }

    Index2D operator+(const Index2D& p) const {
      return Index2D(w+p.w,h+p.h);
    }

    Index2D& operator+=(const Index2D& p) {
      w+=p.w;
      h+=p.h;
      return *this;
    }

    Index2D operator-(const Index2D& p) const {
      return Index2D(w-p.w,h-p.h);
    }

    Index2D& operator-=(const Index2D& p) {
      w-=p.w;
      h-=p.h;
      return *this;
    }
  };
  
  inline Index2D index2d(unsigned int w,unsigned int y) {
    return Index2D(w,y);
  }
}
