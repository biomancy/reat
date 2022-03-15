use std::path::Path;
use std::str::FromStr;

use crate::cli::shared::stranding::Stranding;

pub fn path(rawpath: &str) -> Result<(), String> {
    let path = Path::new(&rawpath);
    if !path.exists() {
        return Err(format!("{} file doesn't exist or there is no permission to read it", rawpath));
    } else {
        Ok(())
    }
}

pub fn writable(rawpath: &str) -> Result<(), String> {
    // TODO: are there any good way to actually check that file is writeable?
    return Ok(());
    // let path = Path::new(&rawpath);
    //
    // if let Some(parent) = path.parent() {
    //     println!("{} {}", parent.display(), parent.exists());
    //     if parent.exists() {
    //         // && !parent.metadata().unwrap().permissions().readonly() {
    //
    //     }
    // }
    // Err(format!("Path {} seems to be not writable", rawpath))
}

pub fn stranding(stranding: &str) -> Result<(), String> {
    match Stranding::from_str(stranding) {
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

#[cfg(test)]
mod test {
    #[test]
    fn stranding() {
        for symbol in ["u", "s", "f", "s/f", "f/s"] {
            assert!(super::stranding(symbol).is_ok());
        }
        for symbol in [".", "r", "uf", "ff", "rr", "+", "-"] {
            assert!(super::stranding(symbol).is_err())
        }
    }

    #[test]
    fn numeric() {
        let validator = super::numeric(10, 12);
        assert!(validator("9").is_err());
        assert!(validator("10").is_ok());
        assert!(validator("12").is_ok());
        assert!(validator("13").is_err());

        let validator = super::numeric(10, 10);
        assert!(validator("9").is_err());
        assert!(validator("10").is_ok());
        assert!(validator("11").is_err());
    }
}
