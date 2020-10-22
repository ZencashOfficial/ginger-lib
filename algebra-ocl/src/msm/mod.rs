mod msm_worker;
pub use msm_worker::*;

mod variable_base;
pub use variable_base::*;

pub fn get_cpu_utilization() -> f64 {
    use std::env;
    use log::error;

    env::var("MSM_CPU_UTILIZATION")
        .and_then(|v| match v.parse() {
            Ok(val) => Ok(val),
            Err(_) => {
                error!("Invalid MSM_CPU_UTILIZATION! Defaulting to 0...");
                Ok(0f64)
            }
        })
        .unwrap_or(0f64)
        .max(0f64)
        .min(1f64)
}
