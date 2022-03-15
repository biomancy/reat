




macro_rules! assort_strands {
    ($items: ident, $func: expr) => {
        $items.unknown.data.retain(|x| match $func(x) {
            Strand::Forward => {
                $items.forward.data.push(x.into());
                false
            }
            Strand::Reverse => {
                $items.reverse.data.push(x.into());
                false
            }
            Strand::Unknown => true,
        });
    };
}

pub(crate) use assort_strands;
