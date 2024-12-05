use wasm_bindgen::prelude::*;

use logicle::{Logicle, LogicleError, LogicleResult};

#[wasm_bindgen]
pub struct LogicleW {
    inner: Logicle,
}

#[wasm_bindgen]
impl LogicleW {
    pub fn validate_parameters(T: f64, W: f64, M: f64, A: f64, bins: i32) -> Result<(), JsValue> {
        Logicle::validate_parameters(T, W, M, A, bins).map_err(|e| JsValue::from(e.to_string()))
    }

    pub fn new(T: f64, W: f64, M: f64, A: f64, bins: i32) -> Result<LogicleW, JsValue> {
        let inner = Logicle::new(T, W, M, A, bins).map_err(|e| JsValue::from(e.to_string()))?;
        Ok(LogicleW { inner })
    }

    pub fn scale(&self, value: f64) -> Result<f64, JsValue> {
        self.inner
            .scale(value)
            .map_err(|e| JsValue::from(e.to_string()))
    }

    pub fn inverse(&self, scale: f64) -> f64 {
        self.inner.inverse(scale)
    }
}

