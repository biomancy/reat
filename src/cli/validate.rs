use std::path::Path;
use std::str::FromStr;

pub fn path(rawpath: &str) -> Result<(), String> {
    let path = Path::new(&rawpath);
    if !path.exists() {
        return Err(format!("{} doesn't exists", rawpath));
    } else {
        Ok(())
    }
}

pub fn writable(rawpath: &str) -> Result<(), String> {
    let path = Path::new(&rawpath);
    if let Some(parent) = path.parent() {
        if parent.exists() {
            return Ok(());
        }
    }
    Err(format!("Path {} seems to be not writable", rawpath))
}

pub fn stranding(stranding: &str) -> Result<(), String> {
    match super::stranding::Stranding::from_str(stranding) {
        Ok(_) => Ok(()),
        Err(x) => Err(x),
    }
}

pub fn numeric<T>(low: T, upper: T) -> impl Fn(&str) -> Result<(), String>
where
    T: FromStr + std::fmt::Display + std::cmp::PartialOrd + Sized,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    move |val: &str| -> Result<(), String> {
        let integer = val.parse::<T>();
        if integer.is_err() {
            return Err(format!("failed to io {}", val));
        };
        let integer = integer.unwrap();

        if integer < low || integer > upper {
            return Err(format!("Value {} is expected to be inside [{}, {}] range", val, low, upper));
        }
        Ok(())
    }
}
