#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray {

public:
    ray() {}
    ray(const point3& origin, const vec3& direction) {
        this->orig = origin;
        this->dir = direction;
    }

    vec3 direction() const { return this->dir; }
    point3 origin() const { return this->orig; }

    point3 at(double t) const {
        return orig + t * dir;
    }

private:
    point3 orig;
    vec3 dir;
};

#endif