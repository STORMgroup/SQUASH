//! The search's token: a value and its span on the read's q-symbol axis, in 8 bytes.

/// A token and its span on the read's q-symbol axis, in 8 bytes: the value and one word holding the start
/// (`32 - TOKEN_SPAN_BITS` = 21 bits, reads up to ~2 million q-symbols) and the span length (`TOKEN_SPAN_BITS` = 11
/// bits). `Token::new` stops the run on a token that does not fit, so no output can differ silently.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Token {
    pub value: u32,
    packed: u32,
}

const TOKEN_SPAN_BITS: u32 = 11;
const TOKEN_START_MAX: i64 = (1 << (32 - TOKEN_SPAN_BITS)) - 1;
const TOKEN_SPAN_MAX: i64 = (1 << TOKEN_SPAN_BITS) - 1;

impl Token {
    /// A token with an i64 span (checked: 0 <= start < 2^21, 0 <= end - start < 2^11).
    pub fn new(value: u32, start: i64, end: i64) -> Token {
        if start < 0 || start > TOKEN_START_MAX || end < start || end - start > TOKEN_SPAN_MAX {
            panic!("token span {}..{} does not fit the 8-byte token (start < 2^{}, end - start < 2^{}); a read this long \
                    is not supported", start, end, 32 - TOKEN_SPAN_BITS, TOKEN_SPAN_BITS);
        }
        Token { value, packed: ((start as u32) << TOKEN_SPAN_BITS) | (end - start) as u32 }
    }
    pub fn s(&self) -> i64 {
        (self.packed >> TOKEN_SPAN_BITS) as i64
    }
    pub fn e(&self) -> i64 {
        self.s() + (self.packed & TOKEN_SPAN_MAX as u32) as i64
    }
    /// The same span with another value.
    pub fn with_value(&self, value: u32) -> Token {
        Token { value, packed: self.packed }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The 8-byte token keeps value, start and end exactly up to its limits, and refuses more.
    #[test]
    fn eight_byte_tokens_round_trip_and_refuse_what_does_not_fit() {
        assert_eq!(std::mem::size_of::<Token>(), 8);
        for &(v, s, e) in &[(0u32, 0i64, 0i64), (7, 0, 9), (u32::MAX, 558_527, 559_416), (3, TOKEN_START_MAX, TOKEN_START_MAX + TOKEN_SPAN_MAX)] {
            let t = Token::new(v, s, e);
            assert_eq!((t.value, t.s(), t.e()), (v, s, e));
            assert_eq!(t.with_value(9).s(), s);
            assert_eq!(t.with_value(9).e(), e);
        }
        assert_ne!(Token::new(1, 5, 9), Token::new(1, 5, 10));
        for &(s, e) in &[(-1i64, 3i64), (TOKEN_START_MAX + 1, TOKEN_START_MAX + 2), (10, 9), (0, TOKEN_SPAN_MAX + 1)] {
            assert!(std::panic::catch_unwind(|| Token::new(1, s, e)).is_err(), "{}..{}", s, e);
        }
    }
}
