;(function () {
  'use strict';
  var MODULE_NAME = 'angular-black-scholes';

  /**
   * JavaScript adopted from Bernt Arne Odegaard's Financial Numerical Recipes
   * http://finance.bi.no/~bernt/gcc_prog/algoritms/algoritms/algoritms.html
   * by Steve Derezinski, CXWeb, Inc.  http://www.cxweb.com
   * Copyright (C) 1998  Steve Derezinski, Bernt Arne Odegaard
   */

  // get dependencies
  // from browserify modules or global variables
  var angular;
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = MODULE_NAME;
    angular = require('angular');
  } else if (typeof window !== 'undefined') {
    angular = window.angular;
  }

  // build module
  angular
    .module(MODULE_NAME, [])
    .service('blackscholes', [
      function () {

        /**
         * Normal distribution function.
         * @param {Number}
         * @returns {Number}
         */
        this.NDF = function (x) {
          return (
            Math.exp(- (x * x) / 2) /
            Math.sqrt(2 * Math.PI)
          );
        };

        /**
         * Cumulative normal distribution function (interpolated, single precision accuracy).
         * http://stackoverflow.com/questions/5259421/cumulative-distribution-function-in-javascript
         * @param {Number}
         * @returns {Number}
         */
        this.CNDFs = function (x) {
          var a1 =  0.31938153;
          var a2 = -0.356563782;
          var a3 =  1.781477937;
          var a4 = -1.821255978;
          var a5 =  1.330274429;
          var k1;
          var k2;
          if (x < 0) return (1 - this.CNDFs(-x));
          else if (x > 6) return 1.0;
          else {
            k1 = 1.0 / (1.0 + 0.2316419 * x);
            k2 = ((((a5 * k1 + a4) * k1 + a3) * k1 + a2) * k1 + a1) * k1;
            return 1.0 - this.NDF(x) * k2;
          }
        };

        /**
         * Cumulative normal distribution function (interpolated, double precision accuracy).
         * http://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
         * @param {Number}
         * @returns {Number}
         */
        this.CNDFd = function (x) {
          var RT2PI = Math.sqrt(2 * Math.PI);
          var SPLIT = 7.07106781186547;
          var N0 = 220.206867912376;
          var N1 = 221.213596169931;
          var N2 = 112.079291497871;
          var N3 = 33.912866078383;
          var N4 = 6.37396220353165;
          var N5 = 0.700383064443688;
          var N6 = 0.0352624965998911;
          var M0 = 440.413735824752;
          var M1 = 793.826512519948;
          var M2 = 637.333633378831;
          var M3 = 296.564248779674;
          var M4 = 86.7807322029461;
          var M5 = 16.064177579207;
          var M6 = 1.75566716318264;
          var M7 = 0.0883883476483184;
          var z = Math.abs(x);
          var answer = 0.0;
          var k1;
          var k2;
          var k3;
          if (z <= 37) {
            k1 = Math.exp(- z * z / 2);
            if (z < SPLIT) {
              k2 = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
              k3 = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
              answer = k1 * k2 / k3;
            } else {
              k3 = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
              answer = k1 / (RT2PI * k3);
            }
          }
          return (x <= 0) ? answer : 1 - answer;
        };

        /**
         * Calculates probability of occuring below and above target price.
         * @param {Number} price
         * @param {Number} target
         * @param {Number} days
         * @param {Number} volatility
         * @returns {Array}
         */
        this.probability = function (price, target, days, volatility) {
          var p = price;
          var q = target;
          var t = days / 365;
          var v = volatility;
          
          var vt = v * Math.sqrt(t);
          var lnpq = Math.log(q / p);
          
          var d1 = lnpq / vt;
          
          var y = Math.floor(1 / (1 + 0.2316419 * Math.abs(d1)) * 100000) / 100000;
          var z = Math.floor(0.3989423 * Math.exp(-(d1 * d1) / 2) * 100000) / 100000;
          
          var y5 = 1.330274 * Math.pow(y, 5);
          var y4 = 1.821256 * Math.pow(y, 4);
          var y3 = 1.781478 * Math.pow(y, 3);
          var y2 = 0.356538 * Math.pow(y, 2);
          var y1 = 0.3193815 * y;
          
          var x = 1 - z * (y5 - y4 + y3 - y2 + y1);
          x = Math.floor(x * 100000) / 100000;
          if (d1 < 0) x = 1 - x;
         
          var pbelow = Math.floor(x * 1000) / 10; 
          var pabove = Math.floor((1 - x) * 1000) / 10;
                     
          return [pbelow, pabove];
        };

        /**
         * Calculates option params by given historical volatility.
         * @param {Boolean} type - option type (true for call, false for put)
         * @param {Number} S - stock price
         * @param {Number} X - strike price
         * @param {Number} rate - risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
         * @param {Number} days - till expiration
         * @param {Number} HV - historical volatility (in percent)
         * @returns {Object} option
         */
        this.getOptionBySigma = function (type, S, X, rate, days, HV) {
          var r = rate / 100;
          var t = days / 365;
          var sigma = HV / 100;

          var sqrt_t = Math.sqrt(t);
          var exp_rt = Math.exp(-r * t);

          var d1 = (Math.log(S / X) + r * t) / (sigma * sqrt_t) + 0.5 * (sigma * sqrt_t);
          var d2 = d1 - sigma * sqrt_t;

          var CNDF_d1 = type ? this.CNDF(d1) : -this.CNDF(-d1);
          var CNDF_d2 = type ? this.CNDF(d2) : -this.CNDF(-d2);
          var NDF_d1 = this.NDF(d1);

          return {
            price: S * CNDF_d1 - X * exp_rt * CNDF_d2,
            delta: CNDF_d1,
            gamma: NDF_d1 / (S * sigma * sqrt_t),
            vega:  S * sqrt_t * NDF_d1,
            theta: ( -(S * sigma * NDF_d1) / (2 * sqrt_t) - r * X * exp_rt * CNDF_d2 ) / 365,
            rho:   X * t * exp_rt * CNDF_d2
          };
        };

        /**
         * Calculates option implied volatility by given price (Newton-Raphson algorithm).
         * http://quant.stackexchange.com/questions/7761/how-to-compute-implied-volatility-calculation
         * @param {Boolean} type - option type (true for call, false for put)
         * @param {Number} S - stock price
         * @param {Number} X - strike price
         * @param {Number} rate - risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
         * @param {Number} days - till expiration
         * @param {Number} price - of the option
         * @returns {Number} iv
         */
        this.getSigmaByPriceNR = function (type, S, X, rate, days, price) {
          var ITERATIONS = 10;
          var ACCURACY = 0.01;

          // initial guess
          var iv = Math.sqrt(2 * Math.PI / (days / 365)) * price / S;

          // iterate
          var bs;
          var dprice;
          for (var i = 0; i < ITERATIONS; i++) {
            bs = this.getOptionBySigma(type, S, X, rate, days, iv * 100);
            dprice = bs.price - price;
            if (Math.abs(dprice) < ACCURACY) return iv;
            else iv -= dprice / bs.vega;
          }

          // failed
          return 0;
        };

        /**
         * Calculates option implied volatility by given price (Bisection algorithm).
         * http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode65.html
         * @param {Boolean} type - option type (true for call, false for put)
         * @param {Number} S - stock price
         * @param {Number} X - strike price
         * @param {Number} rate - risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
         * @param {Number} days - till expiration
         * @param {Number} price - of the option
         * @returns {Number} iv
         */
        this.getSigmaByPriceB = function (type, S, X, rate, days, price) {
          var ACCURACY = 0.001;
          var ivLeft = 0;
          var ivRight = 200;
          var ITERATIONS = Math.round( Math.log( (ivRight - ivLeft) / ACCURACY ) / Math.LN2 + 1 );

          var iv;
          var bs;
          var dprice;

          for (var i = 0; i < ITERATIONS; i++) {
            iv = (ivLeft + ivRight) / 2;
            bs = this.getOptionBySigma(type, S, X, rate, days, iv);
            dprice = bs.price - price;

            if (Math.abs(dprice) < ACCURACY) return iv;
            else if (dprice > 0) ivRight = iv;
            else ivLeft = iv;
          }

          return 0;
        };

        /**
         * Calculates option implied volatility by given price (my own method, including "negative" values).
         * @param {Boolean} type - option type (true for call, false for put)
         * @param {Number} S - stock price
         * @param {Number} X - strike price
         * @param {Number} rate - risk free interest rate (in percent, e.g. 5.25% -> rate=5.25)
         * @param {Number} days - till expiration
         * @param {Number} price - of the option
         * @returns {Number} iv
         */
        this.getSigmaByPrice = function (type, S, X, rate, days, price) {
          var priceIntrinsic = this.getOptionBySigma(type, S, X, rate, days, 0).price;

          if (price === priceIntrinsic) return 0;
          else if (price > priceIntrinsic) return this.getSigmaByPriceB(type, S, X, rate, days, price);
          else return -this.getSigmaByPriceB(type, S, X, rate, days, priceIntrinsic + (priceIntrinsic - price));
        };

        // set defaults
        this.CNDF = this.CNDFs;
        this.iv = this.getSigmaByPrice;
        this.bs = this.getOptionBySigma;
      }
    ]);
})();
