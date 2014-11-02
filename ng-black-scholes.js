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

    .factory('bsFactory', [function () {

      // normal distribution function
      function _NDF (x) {
        return Math.exp(- (x * x) / 2) / Math.sqrt(2 * Math.PI);
      }

      // cumulative normal distribution function
      // interpolated, single precision accuracy
      // http://stackoverflow.com/questions/5259421/cumulative-distribution-function-in-javascript
      function _CNDFs (x) {
        var k1, k2,
            a1 =  0.31938153,
            a2 = -0.356563782,
            a3 =  1.781477937,
            a4 = -1.821255978,
            a5 =  1.330274429;
        
        if (x < 0) return (1 - _CNDFs(-x));
        else if (x > 6) return 1.0;
        else {
          k1 = 1.0 / (1.0 + 0.2316419 * x);
          k2 = ((((a5 * k1 + a4) * k1 + a3) * k1 + a2) * k1 + a1) * k1;
          return 1.0 - _NDF(x) * k2;
        }
      }

      // cumulative normal distribution function
      // interpolated, double precision accuracy
      // http://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
      function _CNDFd (x) {
        var RT2PI = Math.sqrt(2 * Math.PI),
            SPLIT = 7.07106781186547,
            
            N0 = 220.206867912376,
            N1 = 221.213596169931,
            N2 = 112.079291497871,
            N3 = 33.912866078383,
            N4 = 6.37396220353165,
            N5 = 0.700383064443688,
            N6 = 0.0352624965998911,

            M0 = 440.413735824752,
            M1 = 793.826512519948,
            M2 = 637.333633378831,
            M3 = 296.564248779674,
            M4 = 86.7807322029461,
            M5 = 16.064177579207,
            M6 = 1.75566716318264,
            M7 = 0.0883883476483184,

            z = Math.abs(x),
            answer = 0.0,
            k1, k2, k3;

        if (z <= 37) {
          k1 = Math.exp(- z * z / 2);
          
          if (z < SPLIT) {
            k2 = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
            k3 = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
            answer = k1 * k2 / k3;
          }

          else {
            k3 = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
            answer = k1 / (RT2PI * k3);
          }
        }

        return (x <= 0) ? answer : 1 - answer;
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
        HV -- historical volatility (in percent)
      */
      function getOptionBySigma (type, S, X, rate, days, HV) {
        var r = rate / 100,
            t = days / 365,
            sigma = HV / 100,

            sqrt_t = Math.sqrt(t),
            exp_rt = Math.exp(-r * t),

            d1 = (Math.log(S / X) + r * t) / (sigma * sqrt_t) + 0.5 * (sigma * sqrt_t),
            d2 = d1 - sigma * sqrt_t,

            CNDF_d1 = type ? CNDF(d1) : -CNDF(-d1),
            CNDF_d2 = type ? CNDF(d2) : -CNDF(-d2),
            NDF_d1 = NDF(d1);

        return {
          price: S * CNDF_d1 - X * exp_rt * CNDF_d2,
          delta: CNDF_d1,
          gamma: NDF_d1 / (S * sigma * sqrt_t),
          vega:  S * sqrt_t * NDF_d1,
          theta: ( -(S * sigma * NDF_d1) / (2 * sqrt_t) - r * X * exp_rt * CNDF_d2 ) / 365,
          rho:   X * t * exp_rt * CNDF_d2
        };
      }

      // returns option implied volatility by given price
      // Newton-Raphson algorithm
      // http://quant.stackexchange.com/questions/7761/how-to-compute-implied-volatility-calculation
      /*
        INPUTS:
        type -- option type: call -> true, put -> false
        S -- stock price
        X -- strike price
        rate -- risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
        days -- left to expiration
        price -- option price
      */
      function getSigmaByPriceNR (type, S, X, rate, days, price) {
        var ITERATIONS = 10,
            ACCURACY = 0.01,
            iv, i, bs, dprice;

        // initial guess
        iv = Math.sqrt(2 * Math.PI / (days / 365)) * price / S;

        // iterate
        for (i = 0; i < ITERATIONS; i++) {
          bs = getOptionBySigma(type, S, X, rate, days, iv * 100);
          dprice = bs.price - price;
          if (Math.abs(dprice) < ACCURACY) return iv;
          else iv -= dprice / bs.vega;
        }

        // failed
        return 0;
      }

      // returns option implied volatility by given price
      // Bisection algorithm
      // http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode65.html
      /*
        INPUTS:
        type -- option type: call -> true, put -> false
        S -- stock price
        X -- strike price
        rate -- risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
        days -- left to expiration
        price -- option price
      */
      function getSigmaByPriceB (type, S, X, rate, days, price) {
        var ACCURACY = 0.01,
            iv, ivLeft = 0, ivRight = 500,
            ITERATIONS = Math.round( Math.log( (ivRight - ivLeft) / ACCURACY ) / Math.LN2 + 1 ),
            i, bs, dprice;

        for (i = 0; i < ITERATIONS; i++) {
          iv = (ivLeft + ivRight) / 2;
          bs = getOptionBySigma(type, S, X, rate, days, iv);
          dprice = bs.price - price;

          if (Math.abs(dprice) < ACCURACY) return iv;
          else if (dprice > 0) ivRight = iv;
          else ivLeft = iv;
        }

        return 0;
      }

      // returns option implied volatility
      // my own method (including "negative" values)
      /*
        INPUTS:
        type -- option type: call -> true, put -> false
        S -- stock price
        X -- strike price
        rate -- risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
        days -- left to expiration
        price -- option price
      */
      function getSigmaByPrice (type, S, X, rate, days, price) {
        var sigma = getSigmaByPriceB(type, S, X, rate, days, price),
            priceIntrinsic, priceSurrogate;

        if (sigma !== 0) return sigma;
        else {
          priceIntrinsic = getOptionBySigma(type, S, X, rate, days, 0).price;
          priceSurrogate = priceIntrinsic + (priceIntrinsic - price);
          sigma = -getSigmaByPriceB(type, S, X, rate, days, priceSurrogate);
          return sigma;
        }
      }

      var NDF = _NDF, CNDF = _CNDFs;
      return {
        p: probability,
        bs: getOptionBySigma,
        iv: getSigmaByPrice
      };
    }]);

})();
