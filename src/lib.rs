use average::{Estimate, Quantile};
use crate::LogicleError::unknown;

pub type LogicleResult<T> = std::result::Result<T, LogicleError>;

#[derive(Debug)]
pub struct ArgError {
    pub name: String,
    pub description: String,
}

#[derive(Debug)]
pub enum LogicleError {
    unknown(String),
    argument(ArgError),
}


impl LogicleError {
    pub fn new<T>(description: T) -> Self
    where
        T: Into<String>,
    {
        LogicleError::unknown(description.into())
    }

    pub fn arg<A, T>(arg_name: A, description: T) -> Self
    where
        A: Into<String>,
        T: Into<String>,
    {
        LogicleError::argument(ArgError {
            name: arg_name.into(),
            description: description.into(),
        })
    }

    pub fn other<T>(e: T) -> Self
    where
        T: std::error::Error,
    {
        LogicleError::unknown(e.to_string())
    }
}


impl std::fmt::Display for LogicleError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match &self {
            LogicleError::unknown(description) => write!(f, "{}", description),
            LogicleError::argument(arg_error) => write!(f, "{} : {}", arg_error.name, arg_error.description),
        }
    }
}


impl std::error::Error for LogicleError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }

    fn description(&self) -> &str {
        match &self {
            LogicleError::unknown(description) => description,
            LogicleError::argument(error) => &error.description,
        }
    }

    fn cause(&self) -> Option<&dyn std::error::Error> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}


#[derive(Debug, Default, Clone)]
pub struct Logicle {
    params: LogicleParams,
    fast_min_value: f64,
    fast_max_value: f64,
}
#[derive(Debug, Default,Clone)]
struct LogicleParams {
    T: f64,
    W: f64,
    M: f64,
    A: f64,

    a: f64,
    b: f64,
    c: f64,
    d: f64,
    f: f64,

    w: f64,
    x0: f64,
    x1: f64,
    x2: f64,

    xTaylor: f64,
    taylor: Vec<f64>,

    lookup: Vec<f64>,
    bins: i32,

    tolerance: f64,
}

pub struct ExactLogicle {
    params: LogicleParams,
}

impl ExactLogicle {}

#[derive(Debug, Default)]
pub struct Params {
    pub T: f64,
    pub W: f64,
    pub M: f64,
    pub A: f64,
}

const TAYLOR_LENGTH: usize = 1 << 12;
const DEFAULT_TOLERANCE: f64 = 3.0 * f64::EPSILON;


/// Logicle
///
/// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4761345/
/// https://github.com/burtonrj/cytotransform/blob/main/cytotransform/logicle_ext/Logicle.cpp
///
/// Estimate
/// https://github.com/RGLab/flowCore/blob/152b8e4903039aceb32e7393a98ee2dc69f89759/R/AllClasses.R#L4800
impl Logicle {
    pub fn estimate(data: &[f64], M: Option<f64>, A: f64, q: f64) -> LogicleResult<Params> {
        let mut T = f64::MIN;
        for value in data.iter() {
            if *value > T {
                T = *value;
            }
        }
        
        if T <= 0.0 {
            T = 262144.0;
            return Ok(Params { T, W: 0.0, M: 4.5, A });
        }

        let m = match M {
            None => T.log(10.0) + 1.0,
            Some(m) => m,
        };

        let mut W = 0.0;

        let mut estimator = Quantile::new(q);
        for value in data.iter().filter(|&&x| x.is_finite() && x < 0.0) {
            estimator.add(*value);
        }
        let r = estimator.quantile() + f64::EPSILON;
        if !r.is_nan() {
            W = (m - (T / r.abs()).log(10.0)) / 2.0;
            if W < 0.0 {
                return Err(LogicleError::new("W is negative! Try to increase 'M'"));
            }
        }

        if -A > W || A + W > m - W {
            return Err(LogicleError::arg(
                "A",
                format!("Parameter A must be in range [{},{}]", -W, m - 2.0 * W)));
        }

        Ok(Params { T, W, M: m, A })
    }
    pub fn validate_parameters(T: f64, W: f64, M: f64, A: f64, bins: i32) -> Result<(), LogicleError> {
        if T <= 0.0 {
            return Err(LogicleError::arg("T", "T is not positive"));
        }
        if W < 0.0 {
            return Err(LogicleError::arg("W", "W is negative"));
        }
        if M <= 0.0 {
            return Err(LogicleError::arg("M", "M is not positive"));
        }
        if (2.0 * W) > M {
            return Err(LogicleError::arg("W", "W is too large"));
        }
        if -A > W || A + W > M - W {
            return Err(LogicleError::arg("A", "A is too large"));
        }
        if bins > 4096 {
            return Err(LogicleError::arg("bins", "bins max value is 4096"));
        }
        if bins < 8 {
            return Err(LogicleError::arg("bins", "bins min value is 8"));
        }
        Ok(())
    }

    pub fn new(T: f64, W: f64, M: f64, A: f64, bins: i32) -> LogicleResult<Self> {
        Self::validate_parameters(T, W, M, A, bins)?;

        let mut corrected_a = A;
        // if we're going to bin the data make sure that
        // zero is on a bin boundary by adjusting A
        if (bins > 0)
        {
            let mut zero = (W + A) / (M + A);
            zero = (zero * (bins as f64) + 0.5).floor() / (bins as f64);
            corrected_a = (M * zero - W) / (1.0 - zero);
        }

        let mut params = LogicleParams::default();
        params.tolerance = DEFAULT_TOLERANCE;
        params.T = T;
        params.W = W;
        params.M = M;
        params.A = corrected_a;
        params.bins = bins;

        // actual parameters
        // formulas from biexponential paper
        params.w = W / (M + A);
        params.x2 = A / (M + A);
        params.x1 = params.x2 + params.w;
        params.x0 = params.x2 + 2.0 * params.w;
        params.b = (M + A) * std::f64::consts::LN_10;
        params.d = Logicle::solve(params.b, params.w)?;
        let c_a = (params.x0 * (params.b + params.d)).exp();
        let mf_a = (params.b * params.x1).exp() - c_a / (params.d * params.x1).exp();
        params.a = T / ((params.b.exp() - mf_a) - c_a / params.d.exp());
        params.c = c_a * params.a;
        params.f = -mf_a * params.a;

        // use Taylor series near x1, i.e., data zero to
        // avoid round off problems of formal definition
        params.xTaylor = params.x1 + params.w / 4.0;
        // compute coefficients of the Taylor series
        let mut pos_coef = params.a * (params.b * params.x1).exp();
        let mut neg_coef = -params.c / (params.d * params.x1).exp();
        // 16 is enough for full precision of typical scales
        params.taylor = vec![0.0; TAYLOR_LENGTH];
        for i in 0..TAYLOR_LENGTH {
            pos_coef *= params.b / ((i + 1) as f64);
            neg_coef *= -params.d / ((i + 1) as f64);
            params.taylor[i] = pos_coef + neg_coef;
        }
        params.taylor[1] = 0.0; // exact result of Logicle condition

        // params.initialize();
        let mut logicle = Logicle { params, fast_min_value: 0.0, fast_max_value: 0.0 };
        logicle.initialize_bin();
        logicle.initialize_min_max()?;

        Ok(logicle)
    }

    fn initialize_min_max(&mut self) -> LogicleResult<()> {
        self.fast_min_value = self.inverse_fast(0.0)?;
        self.fast_max_value = self.inverse_fast(1.0 - f64::EPSILON)?;
        Ok(())
    }

    // pub fn slope(&self, mut scale: f64) -> f64
    // {
    //     // reflect negative scale regions
    //     if (scale < self.params.x1) {
    //         scale = 2.0 * self.params.x1 - scale;
    //     }
    //
    //     // compute the slope of the biexponential
    //     self.params.a * self.params.b
    //         * (self.params.b * scale).exp() + self.params.c * self.params.d
    //         / (self.params.d * scale).exp()
    // }

    fn solve(b: f64, w: f64) -> LogicleResult<f64> {
        // w == 0 means its really arcsinh
        if w == 0.0 {
            return Ok(b);
        }

        // precision is the same as that of b
        let tolerance = 2.0 * b * f64::EPSILON;

        // based on RTSAFE from Numerical Recipes 1st Edition
        // bracket the root
        let mut d_lo = 0.0;
        let mut d_hi = b;

        // bisection first step
        let mut d = (d_lo + d_hi) / 2.0;
        let mut last_delta = d_hi - d_lo;
        let mut delta;

        // evaluate the f(w,b) = 2 * (ln(d) - ln(b)) + w * (b + d)
        // and its derivative
        let f_b = -2.0 * b.ln() + w * b;
        let mut f = 2.0 * d.ln() + w * d + f_b;
        let mut last_f = f64::NAN;

        for i in 1..20usize {
            // compute the derivative
            let df = 2.0 / d + w;

            // if Newton's method would step outside the bracket
            // or if it isn't converging quickly enough
            if ((d - d_hi) * df - f) * ((d - d_lo) * df - f) >= 0.0
                || (1.9 * f).abs() > (last_delta * df).abs() {
                // take a bisection step
                delta = (d_hi - d_lo) / 2.0;
                d = d_lo + delta;
                if d == d_lo {
                    return Ok(d); // nothing changed, we're done
                }
            } else {
                // otherwise take a Newton's method step
                delta = f / df;
                let t = d;
                d -= delta;
                if d == t {
                    return Ok(d); // nothing changed, we're done
                }
            }
            // if we've reached the desired precision we're done
            if (delta.abs() < tolerance) {
                return Ok(d);
            }

            last_delta = delta;

            // recompute the function
            f = 2.0 * d.ln() + w * d + f_b;
            if f == 0.0 || f == last_f {
                return Ok(d); // found the root or are not going to get any closer
            }

            last_f = f;

            // update the bracketing interval
            if f < 0.0 {
                d_lo = d;
            } else {
                d_hi = d;
            }
        }

        Err(LogicleError::new("exceeded maximum iterations in solve()"))
    }

    fn initialize_bin(&mut self) {
        let mut vec = vec![0.0; (self.params.bins + 1) as usize];
        vec.iter_mut()
            .enumerate()
            .for_each(|(i, v)| *v = self.inverse_exact((i as f64) / (self.params.bins as f64)));
        self.params.lookup = vec;
    }

    pub fn inverse_exact(&self, mut scale: f64) -> f64 {
        // reflect negative scale regions
        let negative = scale < self.params.x1;
        if negative {
            scale = 2.0 * self.params.x1 - scale;
        }

        // compute the biexponential
        let inverse;
        if scale < self.params.xTaylor {
            // near x1, i.e., data zero use the series expansion
            inverse = self.series_biexponential(scale);
        } else {
            // this formulation has better roundoff behavior
            inverse = (self.params.a
                * (self.params.b * scale).exp() + self.params.f)
                - self.params.c / (self.params.d * scale).exp();
        }

        // handle scale for negative values
        if negative {
            -inverse
        } else {
            inverse
        }
    }

    fn series_biexponential(&self, scale: f64) -> f64
    {
        // Taylor series is around x1
        let x = scale - self.params.x1;
        // note that taylor[1] should be identically zero according
        // to the Logicle condition so skip it here
        let mut sum = self.params.taylor[TAYLOR_LENGTH - 1] * x;

        let mut i = TAYLOR_LENGTH - 2;

        while i >= 2 {
            sum = (sum + self.params.taylor[i]) * x;
            i -= 1;
        }

        (sum * x + self.params.taylor[0]) * x
    }

    fn int_scale(&self, value: f64) -> LogicleResult<usize> {
        // binary search for the appropriate bin
        let mut lo = 0;
        let mut hi = self.params.bins;
        while (lo <= hi)
        {
            let mid = (lo + hi) >> 1;
            let key = self.params.lookup[mid as usize];
            if (value < key) {
                hi = mid - 1;
            } else if (value > key) {
                lo = mid + 1;
            } else if (mid < self.params.bins) {
                return Ok(mid as usize);
            } else {
                // equal to table[bins] which is for interpolation only
                return Err(LogicleError::new(format!("IllegalArgument {}", value)));
            }
        }

        // check for out of range
        if (hi < 0 || lo > self.params.bins) {
            return Err(LogicleError::new(format!("IllegalArgument {}", value)));
        }

        Ok((lo - 1) as usize)
    }

    pub fn scale(&self, value: f64) -> LogicleResult<f64> {
        if value <= self.fast_min_value || value >= self.fast_max_value {
            self.scale_exact(value)
        } else {
            self.scale_fast(value)
        }
    }

    pub fn do_scale_fast(&self, value: f64, index: usize) -> f64 {
        // inverse interpolate the table linearly
        let delta = (value - self.params.lookup[index])
            / (self.params.lookup[index + 1] - self.params.lookup[index]);

        (index as f64 + delta) / self.params.bins as f64
    }

    pub fn scale_fast(&self, value: f64) -> LogicleResult<f64> {
        // lookup the nearest value
        let index = self.int_scale(value)?;
        Ok(self.do_scale_fast(value, index))
    }

    pub fn scale_exact(&self, mut value: f64) -> LogicleResult<f64> {
        // handle true zero separately
        if value == 0.0 {
            return Ok(self.params.x1);
        }

        // reflect negative values
        let negative = value < 0.0;
        if negative {
            value = -value;
        }

        // initial guess at solution
        let mut x;
        if value < self.params.f {
            // use linear approximation in the quasi linear region
            x = self.params.x1 + value / self.params.taylor[0];
        } else {
            // otherwise use ordinary logarithm
            x = (value / self.params.a).ln() / self.params.b;
        }

        // try for double precision unless in extended range
        let mut tolerance = self.params.tolerance;
        if x > 1.0 {
            tolerance = self.params.tolerance * x;
        }

        for i in 0..10 {
            // compute the function and its first two derivatives
            let ae2bx = self.params.a * (self.params.b * x).exp();
            let ce2mdx = self.params.c / (self.params.d * x).exp();
            let y;
            if x < self.params.xTaylor {
                // near zero use the Taylor series
                y = self.series_biexponential(x) - value;
            } else {
                // this formulation has better roundoff behavior
                y = (ae2bx + self.params.f) - (ce2mdx + value);
            }
            let abe2bx = self.params.b * ae2bx;
            let cde2mdx = self.params.d * ce2mdx;
            let dy = abe2bx + cde2mdx;
            let ddy = self.params.b * abe2bx - self.params.d * cde2mdx;

            // this is Halley's method with cubic convergence
            let delta = y / (dy * (1.0 - y * ddy / (2.0 * dy * dy)));
            x -= delta;

            // if we've reached the desired precision we're done
            if delta.abs() < tolerance {
                // handle negative arguments
                if (negative) {
                    return Ok(2.0 * self.params.x1 - x);
                } else {
                    return Ok(x);
                }
            }
        }

        Err(LogicleError::new(format!("scale() didn't converge {}", value)))
    }

    pub fn inverse(&self, scale: f64) -> f64 {
        // find the bin
        let x = scale * self.params.bins as f64;
        let index = x.floor() as i32;
        if index < 0 || index >= self.params.bins {
            self.inverse_exact(scale)
        } else {
            self.do_fast_inverse(x, index as usize)
        }
    }

    fn do_fast_inverse(&self, x: f64, index: usize) -> f64 {
        // interpolate the table linearly
        let delta = x - index as f64;

        (1.0 - delta)
            * self.params.lookup[index]
            + delta * self.params.lookup[index + 1]
    }

    pub fn inverse_fast(&self, scale: f64) -> LogicleResult<f64> {
        // find the bin
        let x = scale * self.params.bins as f64;
        let index = x.floor() as i32;
        if index < 0 || index >= self.params.bins {
            return Err(LogicleError::new(format!("IllegalArgument {}", scale)));
        }

        Ok(self.do_fast_inverse(x, index as usize))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_logical() {
        let values = [-10.0, -5.0, -1.0, 0.0, 0.3, 1.0, 3.0, 10.0, 100.0, 999.0];
        let t = 1000.0;
        let w = 1.0;
        let m = 4.0;
        let a = 0.0;

        let expected_scales = [
            0.067574,
            0.147986,
            0.228752,
            0.25,
            0.256384,
            0.271248,
            0.312897,
            0.432426,
            0.739548,
            0.99988997,
        ];

        let logicle = Logicle::new(t, w, m, a, 4096).unwrap();

        let scales = values.iter()
            .map(|value| logicle.scale_fast(*value).unwrap())
            .collect::<Vec<_>>();

        eprintln!("scale_fast {:?}", scales);

        expected_scales.iter().zip(scales.iter()).for_each(|(a, b)| {
            assert!((a - b).abs() < 1e-5);
        });

        let invert_scales = scales.iter()
            .map(|scale| logicle.inverse_fast(*scale).unwrap())
            .collect::<Vec<_>>();

        eprintln!("inverse_fast {:?}", invert_scales);

        invert_scales.iter().zip(values.iter()).for_each(|(a, b)| {
            assert!((a - b).abs() < 1e-5);
        });
    }

    #[test]
    fn test_exact_logical() {
        let values = [-10.0, -5.0, -1.0, 0.0, 0.3, 1.0, 3.0, 10.0, 100.0, 999.0];
        let t = 1000.0;
        let w = 1.0;
        let m = 4.0;
        let a = 0.0;

        let expected_scales = [
            0.067574,
            0.147986,
            0.228752,
            0.25,
            0.256384,
            0.271248,
            0.312897,
            0.432426,
            0.739548,
            0.99988997,
        ];

        let logicle = Logicle::new(t, w, m, a, 4096).unwrap();

        let scales = values.iter()
            .map(|value| logicle.scale_exact(*value).unwrap())
            .collect::<Vec<_>>();

        eprintln!("scale_newton {:?}", scales);

        expected_scales.iter().zip(scales.iter()).for_each(|(a, b)| {
            assert!((a - b).abs() < 1e-5);
        });

        let invert_scales = scales.iter()
            .map(|scale| logicle.inverse_exact(*scale))
            .collect::<Vec<_>>();

        eprintln!("inverse_newton {:?}", invert_scales);

        invert_scales.iter().zip(values.iter()).for_each(|(a, b)| {
            assert!((a - b).abs() < 1e-5);
        });
    }

    #[test]
    fn test_logical() {
        let t = 1000.0;
        let w = 1.0;
        let m = 4.0;
        let a = 0.0;

        let logicle = Logicle::new(t, w, m, a, 4096).unwrap();

        eprintln!("logicle.fast_min_value -- {:?}", logicle.fast_min_value);
        eprintln!("logicle.fast_max_value -- {:?}", logicle.fast_max_value);
 
        assert!(logicle.scale_fast(logicle.fast_min_value - 0.1).is_err());
        let result = logicle.scale_exact(logicle.fast_min_value - 0.1);
        assert!(result.is_ok());
        assert_eq!(logicle.scale_exact(logicle.fast_min_value - 0.1).unwrap(),
                   logicle.scale(logicle.fast_min_value - 0.1).unwrap());

        eprintln!("scale -- fast_min_value - 0.1 {:?}", result.unwrap());

        assert!(logicle.scale_fast(logicle.fast_max_value + 0.1).is_err());
        let result = logicle.scale_exact(logicle.fast_max_value + 0.1);
        assert!(result.is_ok());
        assert_eq!(logicle.scale_exact(logicle.fast_max_value + 0.1).unwrap(),
                   logicle.scale(logicle.fast_max_value + 0.1).unwrap());

        eprintln!("scale -- fast_max_value + 0.1 {:?}", result.unwrap());

        assert!(logicle.inverse_fast(-0.1).is_err());
        let result = logicle.inverse_exact(-0.1);
        
        assert_eq!(logicle.inverse_exact(-0.1), logicle.inverse(-0.1));

        eprintln!("inverse -- result {:?}", result);

        assert!(logicle.inverse_fast(1.1).is_err());
        let result = logicle.inverse_exact(1.1);

        assert_eq!(logicle.inverse_exact(1.1), logicle.inverse(1.1));

        eprintln!("inverse -- result {:?}", result);
    }
}
