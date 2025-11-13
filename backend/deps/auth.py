from __future__ import annotations

import logging
from typing import Optional, Dict, Any

from fastapi import Depends, Header
from sqlalchemy.orm import Session

from db.session import get_db
from db import models


logger = logging.getLogger(__name__)


async def get_current_user(
	db: Session = Depends(get_db),
	x_google_sub: Optional[str] = Header(default=None, alias="X-Google-Sub"),
	x_user_email: Optional[str] = Header(default=None, alias="X-User-Email"),
	x_user_name: Optional[str] = Header(default=None, alias="X-User-Name"),
) -> Dict[str, Any]:
	"""Bind request to a user record.

	- If X-Google-Sub is provided: upsert user by google_sub.
	- If missing: use a single persisted test user (sub='test-sub').
	- Returns lightweight dict with id, google_sub, email, display_name.
	"""

	if x_google_sub:
		user = (
			db.query(models.User)
			.filter(models.User.google_sub == x_google_sub)
			.one_or_none()
		)
		if user is None:
			user = models.User(
				google_sub=x_google_sub,
				email=x_user_email,
				display_name=x_user_name,
			)
			db.add(user)
			db.commit()
			db.refresh(user)
		else:
			updated = False
			if x_user_email and user.email != x_user_email:
				user.email = x_user_email
				updated = True
			if x_user_name and user.display_name != x_user_name:
				user.display_name = x_user_name
				updated = True
			if updated:
				db.commit()

		logger.info(
			"Auth bind: google user sub=%s email=%s name=%s -> id=%s",
			x_google_sub,
			x_user_email,
			x_user_name,
			getattr(user, "id", None),
		)
		return _to_light_user(user)

	# Fallback: single test user persisted once
	test_sub = "test-sub"
	user = (
		db.query(models.User)
		.filter(models.User.google_sub == test_sub)
		.one_or_none()
	)
	if user is None:
		user = models.User(google_sub=test_sub, email="test@example.com", display_name="Test User")
		db.add(user)
		db.commit()
		db.refresh(user)

	logger.info("Auth bind: fallback test user -> id=%s", getattr(user, "id", None))
	return _to_light_user(user)


def _to_light_user(user: models.User) -> Dict[str, Any]:
	return {
		"id": user.id,
		"google_sub": user.google_sub,
		"email": user.email,
		"display_name": user.display_name,
	}


