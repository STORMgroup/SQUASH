//! Shared helpers of the rounds' tests.

use super::*;

pub(super) fn unit(values: &[u32]) -> Vec<Token> {
    values.iter().enumerate().map(|(i, &v)| Token::new(v, i as i64, i as i64 + 1)).collect()
}
