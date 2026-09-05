from __future__ import annotations

import pandas as pd

from assay_context.splits import position_folds


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "variant": [f"A{position}V" for position in range(1, 7) for _ in range(2)],
            "position": [position for position in range(1, 7) for _ in range(2)],
            "primary_endpoint": ["function-impaired", "function-preserved"] * 6,
        }
    )


def test_position_folds_is_group_isolated_and_covers_each_row_once() -> None:
    frame = _frame()
    folds = position_folds(frame, 3)
    test_indices = [index for _, test in folds for index in test]
    assert sorted(test_indices) == list(range(len(frame)))
    assert len(set(test_indices)) == len(frame)
    for train, test in folds:
        assert set(frame.iloc[train]["position"]).isdisjoint(frame.iloc[test]["position"])
