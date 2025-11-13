from __future__ import annotations

import enum
from datetime import datetime
from typing import Optional

import sqlalchemy as sa
from sqlalchemy import ForeignKey, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship
from sqlalchemy.types import TypeDecorator
from sqlalchemy.dialects import postgresql as pg

from .session import Base


class JSONBCompat(TypeDecorator):
	"""JSONB on Postgres, JSON on other dialects (e.g., SQLite)."""
	impl = sa.JSON
	cache_ok = True

	def load_dialect_impl(self, dialect):
		if dialect.name == "postgresql":
			return dialect.type_descriptor(pg.JSONB())
		return dialect.type_descriptor(sa.JSON())


class InferenceStatus(str, enum.Enum):
	pending = "pending"
	processing = "processing"
	completed = "completed"
	failed = "failed"


class ArtifactType(str, enum.Enum):
	predictions_csv = "predictions_csv"
	umap_plot = "umap_plot"
	enriched_h5ad = "enriched_h5ad"
	raw_h5ad = "raw_h5ad"
	log = "log"
	umap_zip = "umap_zip"


class User(Base):
	__tablename__ = "users"

	id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
	google_sub: Mapped[Optional[str]] = mapped_column(sa.String(255), unique=True, nullable=True)
	email: Mapped[Optional[str]] = mapped_column(sa.String(320), nullable=True)
	display_name: Mapped[Optional[str]] = mapped_column(sa.String(255), nullable=True)
	created_at: Mapped[datetime] = mapped_column(sa.DateTime(timezone=True), default=datetime.utcnow, nullable=False)

	inferences: Mapped[list[Inference]] = relationship(back_populates="user")


class Dataset(Base):
	"""Optional datasets catalog, referenced by inferences.dataset_id."""
	__tablename__ = "datasets"

	id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
	name: Mapped[Optional[str]] = mapped_column(sa.String(255), nullable=True)
	description: Mapped[Optional[str]] = mapped_column(sa.Text, nullable=True)
	created_at: Mapped[datetime] = mapped_column(sa.DateTime(timezone=True), default=datetime.utcnow, nullable=False)

	inferences: Mapped[list[Inference]] = relationship(back_populates="dataset")


class Inference(Base):
	__tablename__ = "inferences"

	id: Mapped[str] = mapped_column(sa.String(36), primary_key=True)  # UUID as string to keep codepaths simple
	user_id: Mapped[Optional[int]] = mapped_column(ForeignKey("users.id"), nullable=True, index=True)
	status: Mapped[InferenceStatus] = mapped_column(sa.Enum(InferenceStatus, name="inference_status"), nullable=False, index=True)
	upload_time: Mapped[Optional[datetime]] = mapped_column(sa.DateTime(timezone=True), nullable=True)
	start_time: Mapped[Optional[datetime]] = mapped_column(sa.DateTime(timezone=True), nullable=True)
	end_time: Mapped[Optional[datetime]] = mapped_column(sa.DateTime(timezone=True), nullable=True)
	error_text: Mapped[Optional[str]] = mapped_column(sa.Text, nullable=True)
	model_name: Mapped[Optional[str]] = mapped_column(sa.String(255), nullable=True)
	model_version: Mapped[Optional[str]] = mapped_column(sa.String(255), nullable=True)
	train_args_path: Mapped[Optional[str]] = mapped_column(sa.Text, nullable=True)
	code_commit_sha: Mapped[Optional[str]] = mapped_column(sa.String(64), nullable=True)
	runtime_meta: Mapped[Optional[dict]] = mapped_column(JSONBCompat, nullable=True)
	dataset_id: Mapped[Optional[int]] = mapped_column(ForeignKey("datasets.id"), nullable=True, index=True)

	user: Mapped[Optional[User]] = relationship(back_populates="inferences")
	dataset: Mapped[Optional[Dataset]] = relationship(back_populates="inferences")
	artifacts: Mapped[list[Artifact]] = relationship(back_populates="inference", cascade="all, delete-orphan")


class Artifact(Base):
	__tablename__ = "artifacts"
	__table_args__ = (
		UniqueConstraint("inference_id", "a_type", "storage_path", name="uq_artifact_inference_type_path"),
	)

	id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
	inference_id: Mapped[str] = mapped_column(ForeignKey("inferences.id"), index=True, nullable=False)
	a_type: Mapped[ArtifactType] = mapped_column(sa.Enum(ArtifactType, name="artifact_type"), nullable=False)
	storage_path: Mapped[str] = mapped_column(sa.Text, nullable=False)  # relative path like uploads/{job_id}/...
	mime_type: Mapped[Optional[str]] = mapped_column(sa.String(255), nullable=True)
	size_bytes: Mapped[Optional[int]] = mapped_column(sa.BigInteger, nullable=True)
	sha256: Mapped[Optional[str]] = mapped_column(sa.String(64), nullable=True)
	created_at: Mapped[datetime] = mapped_column(sa.DateTime(timezone=True), default=datetime.utcnow, nullable=False)

	inference: Mapped[Inference] = relationship(back_populates="artifacts")


