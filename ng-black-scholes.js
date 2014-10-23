/* global angular */
;(function () {
  'use strict';

  angular
    .module('black-scholes', [])

    //  JavaScript adopted from Bernt Arne Odegaard's Financial Numerical Recipes
    //  http://finance.bi.no/~bernt/gcc_prog/algoritms/algoritms/algoritms.html
    //  by Steve Derezinski, CXWeb, Inc.  http://www.cxweb.com
    //  Copyright (C) 1998  Steve Derezinski, Bernt Arne Odegaard
    //
    //  This program is free software; you can redistribute it and/or
    //  modify it under the terms of the GNU General Public License
    //  as published by the Free Software Foundation.
     
    //  This program is distributed in the hope that it will be useful,
    //  but WITHOUT ANY WARRANTY; without even the implied warranty of
    //  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    //  GNU General Public License for more details.
    //  http://www.fsf.org/copyleft/gpl.html

    .factory('black-scholes', [function () {

      // normal distribution function
      function NDF (z) {
        return (1.0 / Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * z);
      }

      // cumulative normal distribution function (interpolated)
      // http://stackoverflow.com/questions/5259421/cumulative-distribution-function-in-javascript
      function CNDF (z) {
        var b1, b2, b3, b4, b5, p, c2, a, t, b, n;

        b1 =  0.31938153;
        b2 = -0.356563782;
        b3 =  1.781477937;
        b4 = -1.821255978;
        b5 =  1.330274429;
        p  =  0.2316419;
        c2 =  0.3989423;
        
        a = Math.abs(z);

        if (a > 6.0) return 1.0;

        t = 1.0 / (1.0 + a * p);
        b = c2 * Math.exp(- (z * z) / 2.0);
        n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
        n = 1.0 - b * n;
        if (z < 0.0) n = 1.0 - n;

        return n;
      }  

      // given a decimal number z, returns a string with whole number + fractional string
      // i.e. z = 4.375, returns "4 3/8"
      function fraction (z) {
        var whole = Math.floor(z),
            fract = z - whole,
            thirtytwos = Math.round(fract * 32),
            tmp;

        // if fraction is < 1/64
        if (thirtytwos === 0) return whole + ' ';
        
        // if fraction is > 63/64
        if (thirtytwos === 32) return whole + 1;
       

        // -- 32's non-trivial denominators: 2,4,8,16 -----

        if (thirtytwos / 16 === 1) return whole + ' 1/2';
       
        if (thirtytwos / 8 === 1) return whole + ' 1/4';
        if (thirtytwos / 8 === 3) return whole + ' 3/4';
       
        tmp = thirtytwos / 4;
        if (tmp === Math.floor(tmp)) return whole + ' ' + tmp + '/8';
       
        tmp = thirtytwos / 2;
        if (tmp === Math.floor(tmp)) return whole + ' ' + tmp + '/16';
        

        //
        return whole + ' ' + thirtytwos + '/32';
      }

      // returns probability of occuring below and above target price
      function probability (price, target, days, volatility) {
        var p, q, t, v, vt, lnpq, d1, y, z, y5, y4, y3, y2, y1, x, pbelow, pabove;

        p = price;
        q = target;
        t = days / 365;
        v = volatility;
        
        vt = v * Math.sqrt(t);
        lnpq = Math.log(q / p);
        
        d1 = lnpq / vt;
        
        y = Math.floor(1 / (1 + 0.2316419 * Math.abs(d1)) * 100000) / 100000;
        z = Math.floor(0.3989423 * Math.exp(-(d1 * d1) / 2) * 100000) / 100000;
        
        y5 = 1.330274 * Math.pow(y, 5);
        y4 = 1.821256 * Math.pow(y, 4);
        y3 = 1.781478 * Math.pow(y, 3);
        y2 = 0.356538 * Math.pow(y, 2);
        y1 = 0.3193815 * y;
        
        x = 1 - z * (y5 - y4 + y3 - y2 + y1);
        x = Math.floor(x * 100000) / 100000;
        
        if (d1 < 0) x = 1 - x;
       
        pbelow = Math.floor(x * 1000) / 10; 
        pabove = Math.floor((1 - x) * 1000) / 10;
                   
        return [pbelow, pabove];
      }

      // returns option params by given historical volatility
      /*
        INPUTS:
        type -- option type: call -> true, put -> false
        S -- stock price
        X -- strike price
        rate -- risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
        days -- left to expiration
        HV -- historical volatility
      */
      function option (type, S, X, rate, days, HV) {
        var r = rate / 100,
            t = days / 365,
            sqrt_t   = Math.sqrt(t),
            exp_rt   = Math.exp(-r * t),
            NDF_d1   = NDF(d1),
            CNDF_d2  = type ? CNDF(d2) : -CNDF(-d2),
            delta    = type ? CNDF(d1) : -CNDF(-d1),

            d1 = (Math.log(S / X) + r * t) / (HV * sqrt_t) + 0.5 * (HV * sqrt_t),
            d2 = d1 - (HV * sqrt_t);              

        return {
          price: S * delta - X * exp_rt * CNDF_d2,
          delta: delta,
          gamma: NDF_d1 / (S * HV * sqrt_t),
          vega:  S * sqrt_t * NDF_d1,
          theta: -(S * HV * NDF_d1) / (2 * sqrt_t) - r * X * exp_rt * CNDF_d2,
          rho:   X * t * exp_rt * CNDF_d2
        };
      }

      // returns option implied volatility by given price
      /*
        INPUTS:
        type -- option type: call -> true, put -> false
        S -- stock price
        X -- strike price
        rate -- risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
        days -- left to expiration
        price -- option price
      */
      function iv (type, S, X, rate, days, price) {
        var ITERATIONS = 100,
            ACCURACY = 0.0001,

            r = rate / 100,
            t = days / 365,
            sqrt_t = Math.sqrt(t),

            sigma = (price / S) / (0.398 * sqrt_t),
            i, diff, d1, vega;

        for (i = 0; i < ITERATIONS; i++) {
          diff = price - option(type, S, X, r, sigma, t).price;
          if (Math.abs(diff) < ACCURACY) return sigma;
          else {
            d1 = (Math.log(S / X) + r * t) / (sigma * sqrt_t) + 0.5 * sigma * sqrt_t;
            vega = S * sqrt_t * NDF(d1);
            sigma = sigma + diff / vega;
          }
        }

        // failed to converge
        return -1;
      }

      // factory return
      return {
        fraction: fraction,
        probability: probability,
        option: option,
        iv: iv
      };
    }]);

})();
