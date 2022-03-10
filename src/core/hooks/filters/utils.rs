pub fn maskiter(iter: impl Iterator<Item = bool>) -> Option<Vec<bool>> {
    let mut oneok = false;
    let mask = iter
        .map(|x| {
            oneok = oneok || x;
            x
        })
        .collect();

    match oneok {
        true => Some(mask),
        false => None,
    }
}
